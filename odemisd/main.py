#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on 26 Mar 2012

@author: Éric Piel

Copyright © 2012 Éric Piel, Delmic

This file is part of Open Delmic Microscope Software.

Delmic Acquisition Software is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.

Delmic Acquisition Software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Delmic Acquisition Software. If not, see http://www.gnu.org/licenses/.
'''

from gui import dagui
from odemisd import modelgen
import __version__
import argparse
import logging
import model
import os
import sys

_hwcomponents = set()
def updateMetadata(metadata, parent):
    """
    Update/fill the metadata with all the metadata from all the components
      affecting the given component
    metadata (dict str -> value): metadata
    parent (HwComponent): the component which created the data to which the metadata refers to. 
      Note that the metadata from this very component are not added.
    """
    # find every component which affects the parent
    for comp in _hwcomponents:
        try:
            if parent in comp.affects:
                metadata.update(comp.getMetadata())
        except AttributeError:
            # no affects == empty set
            pass

class BackendContainer(model.Container):
    """
    A normal container which also terminates all the other containers when it
    terminates.
    """
    def __init__(self, name=model.BACKEND_NAME):
        model.Container.__init__(self, name)
        self.sub_containers = set() # to be updated later on
    
    def terminate(self):  
        for container in self.sub_containers:
            try:
                container.terminate()
            except:
                logging.warning("Failed to terminate container %r", container)
                pass
    
        model.Container.terminate(self)
    
    def setMicroscope(self, component):
        self.rootId = component._pyroId

def terminate_all_components(components):
    """
    try to terminate all the components given as much as possible
    components (set of Components): set of components to stop
    """
    for comp in components:
        try:
            comp.terminate()
        except:
            logging.warning("Failed to terminate component '%s'", comp.name)
            pass


BACKEND_RUNNING = "RUNNING"
BACKEND_DEAD = "DEAD"
BACKEND_STOPPED = "STOPPED"
def get_backend_status():
    try:
        microscope = model.getMicroscope()
        if len(microscope.name) > 0:
            return BACKEND_RUNNING
    except:
        if os.path.exists(model.BACKEND_FILE):
            return BACKEND_DEAD
        else:
            return BACKEND_STOPPED
    return BACKEND_DEAD

status_to_xtcode = {BACKEND_RUNNING: 0,
                    BACKEND_DEAD: 1,
                    BACKEND_STOPPED: 2
                    }

# TODO catch kill signal
def run_backend(options):
    model_file = options.model
    try:
        logging.debug("model instantiation file is: %s", model_file.name)
        inst_model = modelgen.get_instantiation_model(model_file)
        logging.info("model has been read successfully")
    except modelgen.ParseError:
        logging.exception("Error while parsing file %s", model_file.name)
        return 127
    
    container = BackendContainer()
    
    try:
        mic, comps, sub_containers = modelgen.instantiate_model(
                                        inst_model, container, 
                                        create_sub_containers=False, # TODO make it possible
                                        dry_run=options.validate)
        # update the model
        _hwcomponents = comps
        container.setMicroscope(mic)
        container.sub_containers |= sub_containers
        logging.info("model has been successfully instantiated")
        logging.debug("model microscope is %s", mic.name) 
        logging.debug("model components are %s", ", ".join([c.name for c in comps])) 
    except:
        logging.exception("When instantiating file %s", model_file.name)
        container.terminate()
        return 127
    
    if options.validate:
        logging.info("model has been successfully validated, exiting")
        terminate_all_components(_hwcomponents)
        container.terminate()
        return 0    # everything went fine
    
    try:
        logging.info("Microscope is now available in container '%s'", model.BACKEND_NAME)
        container.run()
    except:
        # This is coming here in case of signal received when the daemon is running
        logging.exception("When running backend container")
        terminate_all_components(_hwcomponents)
        container.terminate()
        return 127
    
    try:
        terminate_all_components(_hwcomponents)
        container.close()
    except:
        logging.exception("Failed to end the backend container cleanly")
        return 127
    
#    dagui.main(mic)
#    logging.warning("nothing else to do")
    return 0
# This is the cli interface of odemisd, which allows to start the back-end
# It parses the command line and accordingly reads the microscope instantiation
# file, generates a model out of it, and then provides it to the front-end 

def main(args):
    """
    Contains the console handling code for the daemon
    args is the list of arguments passed
    return (int): value to return to the OS as program exit code
    """

    #print args
    # arguments handling 
    parser = argparse.ArgumentParser(description=__version__.name)

    parser.add_argument('--version', action='version', 
                        version=__version__.name + " " + __version__.version + " – " + __version__.copyright)
    dm_grp = parser.add_argument_group('Daemon management')
    dm_grpe = dm_grp.add_mutually_exclusive_group()
    dm_grpe.add_argument("--kill", "-k", dest="kill", action="store_true", default=False,
                        help="Kill the running back-end")
    dm_grpe.add_argument("--check", dest="check", action="store_true", default=False,
                        help="Check for a running back-end (only returns exit code)")
    dm_grpe.add_argument("--daemonize", "-D", action="store_true", dest="daemon",
                         default=False, help="Daemonize the back-end after startup")
    opt_grp = parser.add_argument_group('Options')
    opt_grp.add_argument('--validate', dest="validate", action="store_true", default=False,
                        help="Validate the microscope description file and exit")
    opt_grp.add_argument("--log-level", dest="loglev", metavar="LEVEL", type=int,
                        default=0, help="Set verbosity level (0-2, default = 0)")
    opt_grp.add_argument("--log-target", dest="logtarget", metavar="{auto,stderr,filename}",
                default="auto", help="Specify the log target (auto, stderr, filename)")
    parser.add_argument("model", metavar="file.odm.yaml", nargs='?', type=open, 
                        help="Microscope model instantiation file (*.odm.yaml)")

    options = parser.parse_args(args[1:])
    
    # Set up logging before everything else
    if options.loglev < 0:
        parser.error("log-level must be positive.")
    loglev_names = [logging.WARNING, logging.INFO, logging.DEBUG]
    loglev = loglev_names[min(len(loglev_names) - 1, options.loglev)]
    
    # auto = {odemis.log if daemon, stderr otherwise} 
    if options.logtarget == "auto":
        # default to SysLogHandler ?
        if options.daemon:
            handler = logging.FileHandler("odemis.log")
        else:
            handler = logging.StreamHandler()
    elif options.logtarget == "stderr":
        handler = logging.StreamHandler()
    else:
        handler = logging.FileHandler(options.logtarget)
    logging.getLogger().setLevel(loglev)
    handler.setFormatter(logging.Formatter('%(asctime)s (%(module)s) %(levelname)s: %(message)s'))
    logging.getLogger().addHandler(handler)
    
    # python-daemon is a fancy library but seems to do too many things for us.
    # We just need to contact the backend and see what happens
    
    # Daemon management
    status = get_backend_status()
    if options.kill:
        if status != BACKEND_RUNNING:
            logging.error("No running back-end to kill")
            return 127
        try:
            backend = model.getContainer(model.BACKEND_NAME)
            backend.terminate()
        except:
            logging.error("Failed to stop the back-end")
            return 127
        return 0
    elif options.check:
        logging.info("Status of back-end is %s", status)
        return status_to_xtcode[status]
    
    # check if there is already a backend running
    if status == BACKEND_RUNNING:
        logging.error("Back-end already running, cannot start a new one")
    
    if options.model is None:
        logging.error("No microscope model instantiation file provided")
        return 127
        
    if options.daemon:
        pid = os.fork()
        if pid:
            logging.debug("Daemon started with pid %d", pid)
            return 0
        
    # let's become the backend for real
    return run_backend(options)

if __name__ == '__main__':
    ret = main(sys.argv)
    logging.shutdown() 
    exit(ret)
    
# vim:tabstop=4:shiftwidth=4:expandtab:spelllang=en_gb:spell:
