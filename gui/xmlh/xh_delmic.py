
# This module is used to construct Delmic's custom FoldPanelBar according to
# the definition in the XRC configuration file.
#
# This module is used both by Odemis' GUI and XRCED.

import wx
import wx.xrc as xrc
import odemis.gui.comp.foldpanelbar as fpb
from wx.tools.XRCed.globals import TRACE

class FoldPanelBarXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles
        self.AddStyle('FPB_SINGLE_FOLD', fpb.FPB_SINGLE_FOLD)
        self.AddStyle('FPB_COLLAPSE_TO_BOTTOM', fpb.FPB_COLLAPSE_TO_BOTTOM)
        self.AddStyle('FPB_EXCLUSIVE_FOLD', fpb.FPB_EXCLUSIVE_FOLD)
        self.AddStyle('FPB_HORIZONTAL', fpb.FPB_HORIZONTAL)
        self.AddStyle('FPB_VERTICAL', fpb.FPB_VERTICAL)
        self._isInside = False
        self.current_foldpanelitem = None
        self.spacing = fpb.FPB_DEFAULT_SPACING
        self.left_spacing = fpb.FPB_DEFAULT_LEFTSPACING
        self.right_spacing = fpb.FPB_DEFAULT_RIGHTSPACING

    def CanHandle(self, node):
        # return not self._isInside and self.IsOfClass(node, 'wx.lib.foldpanelbar.FoldPanelBar') or \
        #        self._isInside and self.IsOfClass(node, 'foldpanel')
        return self.IsOfClass(node, 'odemis.gui.comp.foldpanelbar.FoldPanelBar') \
               or self._isInside and self.IsOfClass(node, 'odemis.gui.comp.foldpanelbar.FoldPanelItem') \


    # Process XML parameters and create the object
    def DoCreateResource(self):
        TRACE('DoCreateResource: %s', self.GetClass())

        if self.GetClass() == 'odemis.gui.comp.foldpanelbar.FoldPanelBar':
            w = fpb.FoldPanelBar(self.GetParentAsWindow(),
                                 self.GetID(),
                                 self.GetPosition(),
                                 self.GetSize(),
                                 self.GetStyle(),
                                 self.GetStyle('exstyle'))

            if self.HasParam('spacing'):
                self.spacing = self.GetLong('spacing')

            if self.HasParam('leftspacing'):
                self.left_spacing = self.GetLong('leftspacing')

            if self.HasParam('rightspacing'):
                self.right_spacing = self.GetLong('rightspacing')

            self.SetupWindow(w)
            self._w = w
            old_ins = self._isInside
            self._isInside = True
            # Note: CreateChildren will call this method again
            self.CreateChildren(w, True)
            self._isInside = old_ins

            parent = self._w.GetParent()
            if parent.__class__ == wx.ScrolledWindow:
                parent.EnableScrolling(False, True)
                parent.SetScrollbars(-1, 10, 1, 1)

            return w
        elif self.GetClass() == 'odemis.gui.comp.foldpanelbar.FoldPanelItem':
            item = self._w.AddFoldPanel(self.GetText('label'),
                                        collapsed=self.GetBool('collapsed'),
                                        id=self.GetID())
            self.current_foldpanelitem = item

            n = self.GetParamNode("object")
            wnd = None

            while n:
                #print "Creating Window ", n.GetPropVal('class', "")
                if n.Name != 'object':
                    n = n.Next
                    continue
                wnd = self.CreateResFromNode(n, self.current_foldpanelitem, None)
                if wnd:
                    print self.left_spacing
                    self._w.AddFoldPanelWindow(self.current_foldpanelitem,
                                               wnd,
                                               spacing=self.spacing,
                                               leftSpacing=self.left_spacing,
                                               rightSpacing=self.right_spacing)
                n = n.Next

            # If the last one, was a window ctrl...
            if n and n.Name == 'object' and wnd:
                pass


        wx.CallAfter(self._w.FitBar)

class FoldPanelXmlHandler(xrc.XmlResourceHandler):
    def __init__(self):
        xrc.XmlResourceHandler.__init__(self)
        # Standard styles
        self.AddWindowStyles()
        # Custom styles

    def CanHandle(self, node):
        # return not self._isInside and self.IsOfClass(node, 'wx.lib.foldpanelbar.FoldPanelBar') or \
        #        self._isInside and self.IsOfClass(node, 'foldpanel')
        return self.IsOfClass(node, 'odemis.gui.comp.foldpanelbar.FoldPanelBar') \
               or self.IsOfClass(node, 'odemis.gui.comp.foldpanelbar.FoldPanelItem') \


    # Process XML parameters and create the object
    def DoCreateResource(self):
        pass