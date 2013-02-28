#-*- coding: utf-8 -*-

"""
@author: Rinze de Laat

Copyright © 2013 Rinze de Laat, Delmic

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 2 of the License, or (at your option) any later
version.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.

"""

#===============================================================================
# Test module for exploring and testing the Cairo package
# http://cairographics.org/
#===============================================================================

import unittest
import os
import sys

if os.getcwd().endswith('test'):
    os.chdir('../..')
    print "Working directory changed to", os.getcwd()

print os.getcwd()

import wx
import wx.lib.wxcairo
import cairo

import odemis.gui.test.test_gui
from odemis.gui.xmlh import odemis_get_test_resources

SLEEP_TIME = 100 # Sleep timer in milliseconds
MANUAL = True # If manual is set to True, the window will be kept open at the end
INSPECT = False

def loop():
    app = wx.GetApp()
    if app is None:
        return

    while True:
        wx.CallAfter(app.ExitMainLoop)
        app.MainLoop()
        if not app.Pending():
            break

class TestApp(wx.App):
    def __init__(self):
        odemis.gui.test.test_gui.get_resources = odemis_get_test_resources
        self.test_frame = None
        wx.App.__init__(self, redirect=False)

    def OnInit(self):
        self.test_frame = odemis.gui.test.test_gui.xrccairo_frame(None)
        self.test_frame.SetSize((400, 400))
        self.test_frame.Center()
        self.test_frame.Layout()
        self.test_frame.Show()

        return True

def hex_to_rgb(hex_str):
    hex_str = hex_str[-6:]
    return tuple(int(hex_str[i:i+2], 16) / 255.0 for i in range(0, 6, 2))

def hex_to_rgba(hex_str, af=1.0):
    return hex_to_rgb(hex_str) + (af,)

class CairoPanel(wx.Panel):
    def __init__(self):
        pre = wx.PrePanel()
        # the Create step is done later by XRC.
        self.PostCreate(pre)
        self.Bind(wx.EVT_WINDOW_CREATE, self.OnCreate)

    def OnCreate(self, event):
        self.Unbind(wx.EVT_WINDOW_CREATE)
        # Do all extra initialization here
        self.Bind(wx.EVT_PAINT, self.OnPaint)

    def OnPaint(self, evt):
        dc = wx.BufferedPaintDC(self)
        dc.SetBackground(wx.Brush('#000000'))
        dc.Clear()

        self.go_red(dc)

    def go_red(self, dc):
        print "going red"

        sz = self.GetSize()
        dc.SetPen(wx.Pen("#CCCCCC", 1))

        ctx = wx.lib.wxcairo.ContextFromDC(dc)
        ctx.set_line_width(2)

        color = hex_to_rgba("#8ACD12")
        ctx.set_source_rgba(*color)

        x = y = 0

        while x < sz.width * 2 or y < sz.height * 2:
            x += 20
            y += 20
            ctx.move_to(x, 0)
            ctx.line_to(0, y)

        ctx.stroke()

        image = cairo.ImageSurface.create_from_png("src/odemis/gui/img/example/3-sem.png")
        #w = image.get_width()
        #h = image.get_height()

        ctx.set_source_surface (image, 100, 100)
        ctx.paint ()

        # now draw something with cairo


        #ctx.rel_line_to(-200, 0)
        #ctx.close_path()

        ctx.stroke()

class OwnerDrawnComboBoxTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.app = TestApp()
        loop()

    @classmethod
    def tearDownClass(cls):
        if not MANUAL:
            wx.CallAfter(cls.app.Exit)
        else:
            if INSPECT:
                from wx.lib import inspection
                inspection.InspectionTool().Show()
            cls.app.MainLoop()

    def test_load_image(self):
        pass#self.app.test_frame.cairo_panel.go_red()



if __name__ == "__main__":
    unittest.main()

    #app = TestApp()

    #app.MainLoop()
    #app.Destroy()

