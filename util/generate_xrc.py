#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Embed a XRC file into a Python file.
# python ~/alien/Phoenix/wx/tools/pywxrc.py -p -e -o src/odemis/gui/main_xrc.py.new src/odemis/gui/main.xrc

import glob
import os
import re
import sys
from wx.tools.pywxrc import PythonTemplates, XmlResourceCompiler

PythonTemplates.FILE_HEADER = """\
# This file was automatically generated by pywxrc.
# -*- coding: UTF-8 -*-

import wx
import wx.xrc as xrc

__res = None

def get_resources():
    \"\"\" This function provides access to the XML resources in this module.\"\"\"
    global __res
    if __res is None:
        __init_resources()
    return __res

"""

# Compatible with both wxPython3 and 4
PythonTemplates.CLASS_HEADER = """\
class xrc%(windowName)s(wx.%(windowClass)s):
#!XRCED:begin-block:xrc%(windowName)s.PreCreate
    def PreCreate(self, *args):
        \"\"\" This function is called during the class's initialization.

        Override it for custom setup before the window is created usually to
        set additional window styles using SetWindowStyle() and SetExtraStyle().
        \"\"\"
        pass

#!XRCED:end-block:xrc%(windowName)s.PreCreate

    def __init__(self, parent):
        if wx.MAJOR_VERSION == 3:
            # Two stage creation (see http://wiki.wxpython.org/index.cgi/TwoStageCreation)
            pre = wx.Pre%(windowClass)s()
            self.PreCreate(pre)
            get_resources().LoadOn%(windowClass)s(pre, parent, "%(windowName)s")
            self.PostCreate(pre)
        else:
            wx.%(windowClass)s.__init__(self)
            self.PreCreate()
            get_resources().Load%(windowClass)s(self, parent, "%(windowName)s")

        # Define variables for the controls, bind event handlers
"""

PythonTemplates.INIT_RESOURE_HEADER = """\
# ------------------------ Resource data ----------------------

def __init_resources():
    global __res
    __res = xrc.XmlResource()
"""

PythonTemplates.FILE_AS_STRING = """\
    %(filename)s = b'''\\
%(fileData)s'''

"""

PythonTemplates.XML_AS_STRING = """\
    %(filename)s = u'''\\
%(fileData)s'''

"""

# Fix bytearray
PythonTemplates.ADD_FILE_TO_MEMFS = """\
    wx.MemoryFSHandler.AddFile('XRC/%(memoryPath)s/%(filename)s', bytearray(%(filename)s))
"""

PythonTemplates.ADD_XML_TO_MEMFS = """\
    wx.MemoryFSHandler.AddFile('XRC/%(memoryPath)s/%(filename)s', bytearray(%(filename)s.encode('utf-8')))
"""

class OdemisXmlResourceCompiler(XmlResourceCompiler):

    def NodeContainsFilename(self, node):
        """ Does 'node' contain filename information at all? """

        if node.nodeName == "icon_on":
            return True
            
        return XmlResourceCompiler.NodeContainsFilename(self, node)

    # Fixed version, for label of the menu
    def GenerateWidgetClass(self, windowClass, windowName, topWindow, vars):

        # output the header
        outputList = [self.templates.CLASS_HEADER % locals()]

        # Generate an attribute for each named item in the container
        for widget in topWindow.getElementsByTagName("object"):
            if not self.CheckAssignVar(widget): continue
            widgetClass = widget.getAttribute("class")
            widgetClass = re.sub("^wx", "", widgetClass)
            widgetName = widget.getAttribute("name")
            if widgetName != "" and widgetClass != "":
                vars.append(widgetName)
                if widgetClass == "MenuBar":
                    outputList.append(self.templates.FRAME_MENUBAR_VAR % locals())
                elif widgetClass == "MenuItem":
                    outputList.append(self.templates.FRAME_MENUBAR_MENUITEM_VAR % locals())
                elif widgetClass == "Menu":
                    # Only look directly under for the "label"
                    for e in widget.childNodes:
                        if e.nodeType == e.ELEMENT_NODE and e.tagName == "label":
                            label = e
                            break
                    label = label.childNodes[0].data
                    outputList.append(self.templates.FRAME_MENUBAR_MENU_VAR % locals())
#                 elif widgetClass == "ToolBar":
#                     outputList.append(self.templates.FRAME_TOOLBAR_VAR % locals())
                elif widgetClass == "tool":
                    outputList.append(self.templates.FRAME_TOOLBAR_TOOL_VAR % locals())
                elif widgetClass in ('unknown', 'notebookpage', 'separator', 'sizeritem'):
                    pass
                else:
                    outputList.append(self.templates.CREATE_WIDGET_VAR % locals())

        return outputList

    #-------------------------------------------------------------------

    def GenerateInitResourcesEmbedded(self, resourceFilename, resourceDocument):
        outputList = []
        files = []

        resourcePath = os.path.split(resourceFilename)[0]
        memoryPath = self.GetMemoryFilename(os.path.splitext(os.path.split(resourceFilename)[1])[0])
        resourceFilename = self.GetMemoryFilename(os.path.split(resourceFilename)[1])

        self.ReplaceFilenamesInXRC(resourceDocument.firstChild, files, resourcePath)

        filename = resourceFilename
        fileData = resourceDocument.toxml() # what about this? encoding=resourceDocument.encoding)
        outputList.append(self.templates.XML_AS_STRING % locals())

        for f in files:
            filename = self.GetMemoryFilename(f)
            fileData = self.FileToString(os.path.join(resourcePath, f))
            outputList.append(self.templates.FILE_AS_STRING % locals())

        filename = self.GetMemoryFilename(resourceFilename)
        outputList.append(self.templates.ADD_XML_TO_MEMFS % locals())

        for f in files:
            filename = self.GetMemoryFilename(f)
            outputList.append(self.templates.ADD_FILE_TO_MEMFS % locals())

        outputList.append(self.templates.LOAD_RES_MEMFS % locals())

        return "".join(outputList)


def main(args=None):
    if not args:
        args = sys.argv[1:]

    inputFiles = []
    for arg in args:
        inputFiles += glob.glob(arg)

    embedResources = True
    generateGetText = False
    assignVariables = True
    
    comp = OdemisXmlResourceCompiler()

    try:
        outputFilename = os.path.splitext(args[0])[0] + "_xrc.py"
        comp.MakePythonModule(inputFiles, outputFilename,
                              embedResources, generateGetText,
                              assignVariables)
    except IOError as exc:
        print("%s." % str(exc), file=sys.stderr)
    else:
        print("Resources written to %s." % outputFilename, file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1:])
