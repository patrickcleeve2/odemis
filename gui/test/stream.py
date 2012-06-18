# -*- coding: utf-8 -*-

#===============================================================================
# Test module for Odemis' custom FoldPanelBar in odemis.gui.comp
#===============================================================================

import unittest
import os

if os.getcwd().endswith('test'):
    os.chdir('../..')
    print "Working directory changed to", os.getcwd()

import wx
import odemis.gui.test.test_gui


SLEEP_TIME = 100 # Sleep timer in milliseconds
MANUAL = True # If manual is set to True, the window will be kept open at the end

def odemis_get_resources():
    """ This function provides access to the XML handlers needed for
        non-standard controls defined in the XRC file.
    """
    if odemis.gui.test.test_gui.__res == None:    #pylint: disable=W0212
        from odemis.gui.xmlh.xh_delmic import FoldPanelBarXmlHandler, \
            FixedStreamPanelXmlHandler, CustomStreamPanelXmlHandler
        odemis.gui.test.test_gui.__init_resources() #pylint: disable=W0212
        odemis.gui.test.test_gui.__res.InsertHandler(FoldPanelBarXmlHandler()) #pylint: disable=W0212
        odemis.gui.test.test_gui.__res.InsertHandler(FixedStreamPanelXmlHandler()) #pylint: disable=W0212
        odemis.gui.test.test_gui.__res.InsertHandler(CustomStreamPanelXmlHandler()) #pylint: disable=W0212
    return odemis.gui.test.test_gui.__res #pylint: disable=W0212

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
        odemis.gui.test.test_gui.get_resources = odemis_get_resources
        self.test_frame = None
        wx.App.__init__(self, redirect=False)

    def OnInit(self):
        self.test_frame = odemis.gui.test.test_gui.xrcstream_frame(None)
        self.test_frame.SetSize((400, 400))
        self.test_frame.Center()
        self.test_frame.Layout()
        self.test_frame.Show()
        return True

class FoldPanelBarTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.app = TestApp()
        loop()
        # import wx.lib.inspection
        # wx.lib.inspection.InspectionTool().Show()

    @classmethod
    def tearDownClass(cls):
        if not MANUAL:
            wx.CallAfter(cls.app.Exit)
        else:
            cls.app.MainLoop()

    @classmethod
    def dump_win_tree(cls, window, indent=0):
        if not indent:
            print ""

        for child in window.GetChildren():
            print "."*indent, child.__class__.__name__
            cls.dump_win_tree(child, indent + 2)

    @classmethod
    def has_vertical_scrollbar(cls, window):
        """ Checks if the vertical scroll bar is present by comparing client and
            widget width
        """
        return window.GetClientSize().GetWidth() < window.GetSize().GetWidth()

    @classmethod
    def has_horizontal_scrollbar(cls, window):
        """ Checks if the horizontal scrollbar is present by comparing client and
            widget width
        """
        return window.GetClientSize().GetHeight() < window.GetSize().GetHeight()

    @classmethod
    def build_caption_event(cls, foldpanelitem):

        import odemis.gui.comp.foldpanelbar as ofpb

        cap_bar = foldpanelitem.GetCaptionBar()
        event = ofpb.CaptionBarEvent(ofpb.wxEVT_CAPTIONBAR)
        event.SetEventObject(cap_bar)
        event.SetBar(cap_bar)

        return event

    def test_structure(self):
        wx.MilliSleep(SLEEP_TIME)



if __name__ == "__main__":
    unittest.main()

    #app = TestApp()

    #app.MainLoop()
    #app.Destroy()

