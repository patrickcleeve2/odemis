# -*- coding: utf-8 -*-
'''
Created on 18 Jun 2012

@author: Éric Piel

Copyright © 2012 Éric Piel, Delmic

This file is part of Open Delmic Microscope Software.

Delmic Acquisition Software is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.

Delmic Acquisition Software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with Delmic Acquisition Software. If not, see http://www.gnu.org/licenses/.
'''
from Pyro4.core import oneway
from multiprocessing.process import Process
import Pyro4
import logging
import multiprocessing
import os
import time
import urllib
import weakref

# The special read-only attribute which are duplicated on proxy objects 
class roattribute(property):
    """
    A member of an object which will be cached in the proxy when remotely shared.
    It can be modified only before the object is ever shared. (Technically, it
    can still be written afterwards but the values will not be synchronised
    between the containers).
    """
    # the implementation is just a (python) property with only a different name
    # TODO force to not have setter, but I have no idea how to, override setter?
    pass

def get_roattributes(self):
    """
    list all roattributes of an instance
    """
#    members = inspect.getmembers(self.__class__)
#    return [name for name, obj in members if isinstance(obj, roattribute)]
    klass = self.__class__
    roattributes = []
    for key in dir(klass):
        try:
            if isinstance(getattr(klass, key), roattribute):
                roattributes.append(key)
        except AttributeError:
            continue
       
    return roattributes

def dump_roattributes(self):
    """
    list all the roattributes and their value
    """
    return dict([[name, getattr(self, name)] for name in get_roattributes(self)])

def load_roattributes(self, roattributes):
    """
    duplicate the given roattributes into the instance.
    useful only for a proxy class
    """
    for a, value in roattributes.items():
        setattr(self, a, value)


# Container management functions and class

class ContainerObject(Pyro4.core.DaemonObject):
    """Object which represent the daemon for remote access"""
     
    # it'll never be able to answer back if everything goes fine
    @oneway
    def terminate(self):
        """
        stops the server
        """
        self.daemon.terminate()
    
    def instantiate(self, klass, kwargs):
        """
        instantiate a component and publish it
        klass (class): component class
        kwargs (dict (str -> value)): arguments for the __init__() of the component
        returns the new component instantiated
        """
        return self.daemon.instantiate(klass, kwargs)
    
    def getObject(self, objectId):
        """
        return a registered object from its object id.
        This is mostly useful to get a proxy with all the information about the
        object automatically.
        It will raise an exception if the object is not registered.
        """
        assert isinstance(objectId, basestring)
        return self.daemon.objectsById[objectId]
    
    def getRoot(self):
        """
        returns the root object, if it has been defined in the container
        """
        return self.getObject(self.daemon.rootId)
    
# Basically a wrapper around the Pyro Daemon 
class Container(Pyro4.core.Daemon):
    def __init__(self, name):
        """
        name: name of the container (must be unique)
        """
        assert not "/" in name
        # TODO put all the sockets in some specified directory on the fs /var/odemis/?
        self.ipc_name = urllib.quote(name) + ".ipc"
        
        if os.path.exists(self.ipc_name):
            try:
                os.remove(self.ipc_name)
                logging.warning("The file '%s' was deleted to create container '%s'.", self.ipc_name, name)
            except OSError:
                logging.error("Impossible to delete file '%s', needed to create container '%s'.", self.ipc_name, name)

        Pyro4.Daemon.__init__(self, unixsocket=self.ipc_name, interface=ContainerObject)
        
        # To be set by the user of the container 
        self.rootId = None # objectId of a "Root" component
        
    def run(self, *args, **kwargs):
        """
        runs and serve the objects registered in the container.
        returns only when .terminate() is called
        """
        # wrapper to requestLoop() just because the name is strange
        self.requestLoop(*args, **kwargs)

    def terminate(self):
        """
        stops the server
        """
        # wrapper to shutdown(), in order to be more consistent with the vocabulary
        self.shutdown()
        # All the cleaning is done in the original thread, after the run()
    
    def close(self):
        """
        Cleans up everything behind, once the container is already done with running
        """
        # unregister every object still around, to be sure everything gets
        # deallocated from the memory (but normally, it's up to the client to
        # terminate() every component before)
        for obj in self.objectsById.values():
            if hasattr(obj, "_unregister"):
                try:
                    obj._unregister()
                except Exception:
                    logging.exception("Failed to unregister object %s when terminating container", str(obj))
            else:
                self.unregister(obj)
        
        Pyro4.Daemon.close(self)

    def instantiate(self, klass, kwargs):
        """
        instantiate a Component and publish it
        klass (class): component class
        kwargs (dict (str -> value)): arguments for the __init__() of the component
        returns the new component instantiated
        """
        kwargs["daemon"] = self # the component will auto-register
        comp = klass(**kwargs)
        return comp

# helper functions
def getContainer(name):
    """
    returns (a proxy to) the container with the given name
    raises an exception if no such container exist
    """
    # the container is the default pyro daemon at the address named by the container
    container = Pyro4.Proxy("PYRO:Pyro.Daemon@./u:"+urllib.quote(name)+".ipc",
                            oneways=["terminate"])
    # A proxy doesn't connect until the first remote call, check the connection
    container.ping() # raise an exception if connection fails
    return container

def getObject(container_name, object_name):
    """
    returns (a proxy to) the object with the given name in the given container
    raises an exception if no such object or container exist
    """
    container = getContainer(container_name)
    return container.getObject(urllib.quote(object_name))

def createNewContainer(name):
    """
    creates a new container in an independent and isolated process
    returns the (proxy to the) new container
    """
    # create a container separately
    isready = multiprocessing.Event()
    p = Process(name="Container "+name, target=_manageContainer, args=(name,isready))
    p.start()
    if not isready.wait(3): # wait maximum 3s
        logging.error("Container %s is taking too long to get ready", name)
        raise IOError("Container creation timeout")

    # connect to the new container
    return getContainer(name)
 
def createInNewContainer(container_name, klass, kwargs):
    """
    creates a new component in a new container
    container_name (string)
    klass (class): component class
    kwargs (dict (str -> value)): arguments for the __init__() of the component
    returns the (proxy to the) new component
    """
    container = createNewContainer(container_name)
    return container.instantiate(klass, kwargs)
 
def _manageContainer(name, isready=None):
    """
    manages the whole life of a container, from birth till death
    name (string)
    isready (Event): set when the container is (almost) ready to publish objects
    """
    container = Container(name)
    if isready is not None:
        isready.set()
    container.run()
    container.close()
    
# Special functions and class to manage method/function with weakref
# wxpython.pubsub has something similar 

class WeakRefLostError(Exception):
    pass

class WeakMethodBound(object):
    def __init__(self, f):
        self.f = f.im_func
        self.c = weakref.ref(f.im_self)
        # cache the hash so that it's the same after deref'd
        self.hash = hash(f.im_func) + hash(f.im_self)
        
    def __call__(self, *arg, **kwargs):
        ins = self.c() 
        if ins == None:
            raise WeakRefLostError, 'Method called on dead object'
        return self.f(ins, *arg, **kwargs)
        
    def __hash__(self):
        return self.hash
    
    def __eq__(self, other):
        try:
            return (type(self) is type(other) and self.f == other.f
                    and self.c() == other.c())
        except:
            return False
        
#    def __ne__(self, other):
#        return not self == other

class WeakMethodFree(object):
    def __init__(self, f):
        self.f = weakref.ref(f)
        # cache the hash so that it's the same after deref'd
        self.hash = hash(f)
        
    def __call__(self, *arg, **kwargs):
        fun = self.f()
        if fun == None:
            raise WeakRefLostError, 'Function no longer exist'
        return fun(*arg, **kwargs)
        
    def __hash__(self):
        return self.hash
    
    def __eq__(self, other):
        try:
            return type(self) is type(other) and self.f() == other.f()
        except:
            return False
        
#    def __ne__(self, other):
#        return not self == other

def WeakMethod(f):
    try:
        f.im_func
    except AttributeError:
        return WeakMethodFree(f)
    return WeakMethodBound(f)
