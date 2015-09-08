# -*- coding: utf-8 -*-
"""
Created on 8 Feb 2012

:author: Éric Piel
:copyright: © 2012 Éric Piel, Delmic

.. license::

    This file is part of Odemis.

    Odemis is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License version 2 as published by the Free
    Software Foundation.

    Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along with
    Odemis. If not, see http://www.gnu.org/licenses/.

"""

from __future__ import division

from abc import abstractmethod, ABCMeta
import logging
from odemis import gui, model
from odemis.acq.stream import OpticalStream, EMStream, SpectrumStream, StaticStream
from odemis.gui import BG_COLOUR_LEGEND, FG_COLOUR_LEGEND
from odemis.gui.comp import miccanvas
from odemis.gui.comp.canvas import CAN_DRAG, CAN_FOCUS
from odemis.gui.comp.legend import InfoLegend, AxisLegend
from odemis.gui.img.data import getico_blending_goalBitmap
from odemis.gui.model import CHAMBER_VACUUM, CHAMBER_UNKNOWN
from odemis.gui.util import call_in_wx_main
from odemis.gui.util.raster import rasterize_line
from odemis.model import NotApplicableError
from odemis.util import units, spectrum
import random
import wx


class ViewPort(wx.Panel):

    # Default classes for the canvas and the legend. These may be overridden
    # in subclasses
    canvas_class = None
    bottom_legend_class = None
    left_legend_class = None

    def __init__(self, *args, **kwargs):
        """Note: The MicroscopeViewport is not fully initialised until setView()
        has been called.
        """

        hide_bottom_legend = kwargs.pop('hide_bottom_legend', False)
        hide_left_legend = kwargs.pop('hide_left_legend', False)

        wx.Panel.__init__(self, *args, **kwargs)

        self._microscope_view = None  # model.MicroscopeView
        self._tab_data_model = None  # model.MicroscopyGUIData

        # Keep track of this panel's pseudo focus
        self._has_focus = False

        font = wx.Font(8, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL)
        self.SetFont(font)
        self.SetBackgroundColour(BG_COLOUR_LEGEND)
        self.SetForegroundColour(FG_COLOUR_LEGEND)

        # This attribute can be used to track the (GridBag) sizer position of the viewport (if any)
        self.sizer_pos = None

        # main widget
        self.canvas = self.canvas_class(self)

        # Put all together (canvas + legend)
        self.bottom_legend = None
        self.left_legend = None

        main_sizer = wx.BoxSizer(wx.VERTICAL)
        if (
                self.bottom_legend_class and not hide_bottom_legend and
                self.left_legend_class and not hide_left_legend
        ):
            self.bottom_legend = self.bottom_legend_class(self)
            self.left_legend = self.left_legend_class(self, orientation=wx.VERTICAL)

            grid_sizer = wx.GridBagSizer()
            grid_sizer.Add(self.canvas, pos=(0, 1), flag=wx.EXPAND)
            grid_sizer.Add(self.bottom_legend, pos=(1, 1), flag=wx.EXPAND)
            grid_sizer.Add(self.left_legend, pos=(0, 0), flag=wx.EXPAND)

            filler = wx.Panel(self)
            filler.SetBackgroundColour(BG_COLOUR_LEGEND)
            grid_sizer.Add(filler, pos=(1, 0), flag=wx.EXPAND)

            grid_sizer.AddGrowableRow(0, 1)
            grid_sizer.AddGrowableCol(1, 1)
            # grid_sizer.RemoveGrowableCol(0)

            # Focus the view when a child element is clicked
            self.bottom_legend.Bind(wx.EVT_LEFT_DOWN, self.OnChildFocus)
            self.left_legend.Bind(wx.EVT_LEFT_DOWN, self.OnChildFocus)

            main_sizer.Add(grid_sizer, 1, border=2, flag=wx.EXPAND | wx.ALL)
        elif self.bottom_legend_class and not hide_bottom_legend:
            main_sizer.Add(self.canvas, proportion=1, border=2,
                           flag=wx.EXPAND | wx.TOP | wx.LEFT | wx.RIGHT)
            # It's made of multiple controls positioned via sizers
            # TODO: allow the user to pick which information is displayed
            # in the legend
            self.bottom_legend = self.bottom_legend_class(self)
            self.bottom_legend.Bind(wx.EVT_LEFT_DOWN, self.OnChildFocus)

            main_sizer.Add(self.bottom_legend, proportion=0, border=2,
                           flag=wx.EXPAND | wx.BOTTOM | wx.LEFT | wx.RIGHT)
        elif self.left_legend_class:
            raise NotImplementedError("Only left legend not handled")
        else:
            main_sizer.Add(self.canvas, 1, border=2, flag=wx.EXPAND | wx.ALL)

        self.SetSizerAndFit(main_sizer)
        main_sizer.Fit(self)
        self.SetAutoLayout(True)

        self.Bind(wx.EVT_CHILD_FOCUS, self.OnChildFocus)
        self.Bind(wx.EVT_SIZE, self.OnSize)

    def __str__(self):
        return "{0} {2} {1}".format(
            self.__class__.__name__,
            self._microscope_view.name.value if self._microscope_view else "",
            id(self))

    __repr__ = __str__

    @property
    def microscope_view(self):
        return self._microscope_view

    def clear(self):
        self.canvas.clear()
        if self.bottom_legend:
            self.bottom_legend.clear()
        if self.left_legend:
            self.left_legend.clear()
        self.Refresh()

    def setView(self, microscope_view, tab_data):
        raise NotImplementedError

    ################################################
    # Panel control
    ################################################

    def ShowLegend(self, show):
        """ Show or hide the merge slider """
        self.bottom_legend.Show(show)

    def HasFocus(self, *args, **kwargs):
        return self._has_focus is True

    def SetFocus(self, focus):
        """ Set the focus on the viewport according to the focus parameter.

        focus:  A boolean value.

        """

        logging.debug(["Removing focus from %s", "Setting focus to %s"][focus], id(self))

        self._has_focus = focus
        if focus:
            self.SetBackgroundColour(gui.BORDER_COLOUR_FOCUS)
        else:
            self.SetBackgroundColour(gui.BORDER_COLOUR_UNFOCUS)

    ################################################
    # GUI Event handling
    ################################################

    def OnChildFocus(self, evt):
        """ Give the focus to the view if one of the child widgets is clicked """
        if self._microscope_view and self._tab_data_model:
            # This will take care of doing everything necessary
            # Remember, the notify method of the vigilant attribute will
            # only fire if the values changes.
            self._tab_data_model.focussedView.value = self._microscope_view

        evt.Skip()

    def OnSize(self, evt):
        evt.Skip()  # processed also by the parent

    def Disable(self, *args, **kwargs):
        logging.debug("Disabling %s", self.canvas)
        wx.Panel.Disable(self, *args, **kwargs)
        self.canvas.Disable(*args, **kwargs)

    def Enable(self, *args, **kwargs):
        logging.debug("Enabling %s", self.canvas)
        wx.Panel.Enable(self, *args, **kwargs)
        self.canvas.Enable(*args, **kwargs)


class CameraViewport(ViewPort):
    """ Simple viewport for displaying a video feed, with any added parameters """
    canvas_class = miccanvas.BitmapCanvas


class MicroscopeViewport(ViewPort):
    """ A panel that shows a microscope view and its legend(s)

    This is a generic class, that should be inherited by more specific classes.

    """
    canvas_class = miccanvas.DblMicroscopeCanvas
    bottom_legend_class = InfoLegend

    def __init__(self, *args, **kwargs):
        """Note: The MicroscopeViewport is not fully initialised until setView()
        has been called.
        """
        # Call parent constructor at the end, because it needs the legend panel
        super(MicroscopeViewport, self).__init__(*args, **kwargs)

        if self.bottom_legend:
            # Bind on EVT_SLIDER to update even while the user is moving
            self.bottom_legend.Bind(wx.EVT_LEFT_UP, self.OnSlider)
            self.bottom_legend.Bind(wx.EVT_SLIDER, self.OnSlider)

        # Find out screen pixel density for "magnification" value. This works
        # only with monitors/OS which report correct values. It's unlikely to
        # work with a projector!
        # The 24" @ 1920x1200 screens from Dell have an mpp value of 0.00027 m/px
        self._mpp_screen = 1e-3 * wx.DisplaySizeMM()[0] / wx.DisplaySize()[0]

        # This attribute is set to True if this object (i.e. 'self') was responsible for changing
        # the HFW value
        self.self_set_hfw = False
        self._fov_va = None  # (hardware) VA to follow for HFW

    def setView(self, microscope_view, tab_data):
        """
        Set the microscope view that this viewport is displaying/representing
        *Important*: Should be called only once, at initialisation.

        :param microscope_view:(model.MicroscopeView)
        :param tab_data: (model.MicroscopyGUIData)
        """

        # This is a kind of a kludge, as it'd be best to have the viewport
        # created after the microscope view, but they are created independently
        # via XRC.
        assert(self._microscope_view is None)

        # import traceback
        # traceback.print_stack()

        self._microscope_view = microscope_view
        self._tab_data_model = tab_data

        # TODO: Center to current view position, with current mpp
        microscope_view.mpp.subscribe(self._onMPP, init=True)

        # set/subscribe merge ratio
        microscope_view.merge_ratio.subscribe(self._onMergeRatio, init=True)

        # subscribe to image, to update legend on stream tree/image change
        microscope_view.lastUpdate.subscribe(self._onImageUpdate, init=True)

        # By default, cannot focus, unless the child class allows it
        self.canvas.abilities.discard(CAN_FOCUS)

        # canvas handles also directly some of the view properties
        self.canvas.setView(microscope_view, tab_data)

    ################################################
    #  Panel control
    ################################################

    def ShowMergeSlider(self, show):
        """ Show or hide the merge slider """
        if self.bottom_legend:
            self.bottom_legend.bmp_slider_left.Show(show)
            self.bottom_legend.merge_slider.Show(show)
            self.bottom_legend.bmp_slider_right.Show(show)

    def UpdateHFWLabel(self):
        """ Physical width of the display"""
        if not self._microscope_view or not self.bottom_legend:
            return
        hfw = self._microscope_view.mpp.value * self.GetClientSize()[0]
        hfw = units.round_significant(hfw, 4)
        label = u"HFW: %s" % units.readable_str(hfw, "m", sig=3)
        self.bottom_legend.set_hfw_label(label)

    def UpdateMagnification(self):
        # Total magnification
        mag = self._mpp_screen / self._microscope_view.mpp.value
        label = u"Mag: × %s" % units.readable_str(units.round_significant(mag, 3))

        # Gather all different image mpp values
        mpps = set()
        for im in self._microscope_view.stream_tree.getImages():
            try:
                mpps.add(im.metadata[model.MD_PIXEL_SIZE][0])
            except KeyError:
                pass

        # If there's only one mpp value (i.e. there's only one image, or they
        # all have the same mpp value), indicate the digital zoom.
        if len(mpps) == 1:
            mpp_im = mpps.pop()
            # mag_im = self._mpp_screen / mpp_im  # as if 1 im.px == 1 sc.px
            mag_dig = mpp_im / self._microscope_view.mpp.value
            label += u" (Digital: × %s)" % units.readable_str(units.round_significant(mag_dig, 2))

        if self.bottom_legend:
            self.bottom_legend.set_mag_label(label)

    ################################################
    #  VA handling
    ################################################

    @call_in_wx_main
    def _onMergeRatio(self, val):
        # round is important because int can cause unstable value
        # int(0.58*100) = 57
        if self.bottom_legend:
            self.bottom_legend.merge_slider.SetValue(round(val * 100))

    @call_in_wx_main
    def _onMPP(self, mpp):
        if self.bottom_legend:
            self.bottom_legend.scale_win.SetMPP(mpp)
            self.UpdateHFWLabel()
            self.UpdateMagnification()
            # the MicroscopeView will send an event that the view has to be redrawn

    def _checkMergeSliderDisplay(self):
        """
        Update the MergeSlider display and icons depending on the state
        """

        if not self.bottom_legend:
            return

        # MergeSlider is displayed if:
        # * Root operator of StreamTree accepts merge argument
        # * (and) Root operator of StreamTree has >= 2 images
        if (
                "merge" in self._microscope_view.stream_tree.kwargs
                and len(self._microscope_view.stream_tree) >= 2
        ):
            streams = self._microscope_view.getStreams()
            all_opt = all(isinstance(s, OpticalStream) for s in streams)

            # If all images are optical, assume they are merged using screen blending and no
            # merge ratio is required
            if all_opt:
                self.ShowMergeSlider(False)
            else:
                # TODO: How is the order guaranteed? (Left vs Right)
                # => it should be done in the MicroscopeView when adding a stream
                # For now, special hack for the MicroscopeCanvas which always sets
                # the EM image as "right" (ie, it's drawn last).
                # If there is SEM and Spectrum, the spectrum image is always
                # set as "right" (ie, it's drawn last).
                if (
                        any(isinstance(s, EMStream) for s in streams)
                        and any(isinstance(s, OpticalStream) for s in streams)
                ):
                    self.bottom_legend.set_stream_type(wx.LEFT, OpticalStream)
                    self.bottom_legend.set_stream_type(wx.RIGHT, EMStream)
                elif (
                        any(isinstance(s, EMStream) for s in streams)
                        and any(isinstance(s, SpectrumStream) for s in streams)
                ):
                    self.bottom_legend.set_stream_type(wx.LEFT, EMStream)
                    self.bottom_legend.set_stream_type(wx.RIGHT, SpectrumStream)
                else:
                    sc = self._microscope_view.stream_tree[0]
                    self.bottom_legend.set_stream_type(wx.LEFT, sc.__class__)

                    sc = self._microscope_view.stream_tree[1]
                    self.bottom_legend.set_stream_type(wx.RIGHT, sc.__class__)

                self.ShowMergeSlider(True)
        else:
            self.ShowMergeSlider(False)

    @call_in_wx_main
    def _onImageUpdate(self, _):
        self._checkMergeSliderDisplay()

        # magnification might have changed (eg, image with different binning)
        self.UpdateMagnification()

    ################################################
    # GUI Event handling
    ################################################

    def OnSlider(self, evt):
        """
        Merge ratio slider
        """
        if self._microscope_view is None or not self.bottom_legend:
            return

        val = self.bottom_legend.merge_slider.GetValue() / 100
        self._microscope_view.merge_ratio.value = val
        evt.Skip()

    def OnSize(self, evt):
        evt.Skip()  # processed also by the parent
        self.UpdateHFWLabel()

    def OnSliderIconClick(self, evt):
        evt.Skip()

        if self._microscope_view is None or not self.bottom_legend:
            return

        if evt.GetEventObject() == self.bottom_legend.bmp_slider_left:
            self.bottom_legend.merge_slider.set_to_min_val()
        else:
            self.bottom_legend.merge_slider.set_to_max_val()

        val = self.bottom_legend.merge_slider.GetValue() / 100
        self._microscope_view.merge_ratio.value = val
        evt.Skip()

    # END Event handling

    def track_view_hfw(self, fov_va):
        """ Link the field of view (width) of the view with the field of view of the hardware

        :param fov_va: (FloatVA)

        """

        logging.info("Tracking mpp on %s" % self)
        self._fov_va = fov_va
        # The view FoV changes either when the mpp changes or on resize,
        # but resize typically causes an update of the mpp (to keep the FoV)
        # so no need to listen to resize.
        self.microscope_view.mpp.subscribe(self._on_mpp_set_hfw)
        fov_va.subscribe(self._on_hfw_set_mpp, init=True)

    def _on_hfw_set_mpp(self, hfw):
        """ Change the mpp value of the MicroscopeView when the HFW changes

        We set the mpp value of the MicroscopeView by assigning the microscope's hfw value to
        the Canvas' hfw value, which will cause the the Canvas to calculate a new mpp value
        and assign it to View's mpp attribute.

        """

        # If this ViewPort was not responsible for updating the hardware HFW, update the MPP value
        # of the canvas (which is done in the `horizontal_field_width` setter)
        if not self.self_set_hfw:
            logging.info("Calculating mpp from hfw for viewport %s" % self)
            self.canvas.horizontal_field_width = hfw

        self.self_set_hfw = False

    def _on_mpp_set_hfw(self, _):
        """ Set the microscope's hfw when the MicroscopeView's mpp value changes

        The canvas calculates the new hfw value.
        """

        logging.info("Calculating hfw from mpp for viewport %s" % self)
        hfw = self.canvas.horizontal_field_width

        try:
            # TODO: Test with a simulated SEM that has HFW choices
            choices = self._fov_va.choices
            # Get the choice that matches hfw most closely
            hfw = min(choices, key=lambda choice: abs(choice - hfw))
        except NotApplicableError:
            hfw = self._fov_va.clip(hfw)

        # Indicate that this object was responsible for updating the hardware's HFW, so it won't
        # get updated again in `_on_hfw_set_mpp`
        self.self_set_hfw = True
        self._fov_va.value = hfw


class OverviewViewport(MicroscopeViewport):
    """ A Viewport containing a downscaled overview image of the loaded sample """

    canvas_class = miccanvas.OverviewCanvas
    bottom_legend_class = InfoLegend

    def __init__(self, *args, **kwargs):
        super(OverviewViewport, self).__init__(*args, **kwargs)
        self.Parent.Bind(wx.EVT_SIZE, self.OnSize)

    def OnSize(self, evt):
        # TODO: this can be avoided by just setting a different minimum mpp

        if self.canvas.horizontal_field_width < 10e-3:
            self.canvas.horizontal_field_width = 10e-3
            logging.debug("Canvas HFW too small! Setting it to %s", 10e-3)

        super(OverviewViewport, self).OnSize(evt)
        self.canvas.fit_view_to_content(True)

    def setView(self, microscope_view, tab_data):
        """ Attach the MicroscopeView associated with the overview """

        super(OverviewViewport, self).setView(microscope_view, tab_data)

        self.canvas.point_select_overlay.p_pos.subscribe(self._on_position_select)
        # Only allow moving when chamber is under vacuum
        tab_data.main.chamberState.subscribe(self._on_chamber_state_change, init=True)

    def _on_chamber_state_change(self, chamber_state):
        """ Watch position changes in the PointSelectOverlay if the chamber is ready """

        # If state is unknown, it's probably going to be unknown forever, so
        # we have to allow (and in the worst case the user will be able to move
        # while the chamber is opened)
        if (chamber_state in {CHAMBER_VACUUM, CHAMBER_UNKNOWN} and
                self._microscope_view.has_stage()):
            self.canvas.point_select_overlay.activate()
        else:
            self.canvas.point_select_overlay.deactivate()

    def _on_position_select(self, p_pos):
        """ Set the physical view position
        """
        if self._tab_data_model:
            if self._microscope_view.has_stage():
                self._microscope_view.moveStageTo(p_pos)


class LiveViewport(MicroscopeViewport):
    """
    Used to display live streams on Secom and Delphi.
    The main difference is the handling of the pause state, which prevents
    stage move and indicate it via an icon.
    """

    def __init__(self, *args, **kwargs):
        super(LiveViewport, self).__init__(*args, **kwargs)
        self._orig_abilities = set()

    def setView(self, microscope_view, tab_data):
        # Must be before calling the super, as the super drops CAN_FOCUS automatically
        self._orig_abilities = self.canvas.abilities & {CAN_DRAG, CAN_FOCUS}
        super(LiveViewport, self).setView(microscope_view, tab_data)
        tab_data.streams.subscribe(self._on_stream_change)
        microscope_view.stream_tree.should_update.subscribe(self._on_stream_play,
                                                            init=True)

    def _on_stream_play(self, is_playing):
        """
        Called whenever view contains a stream playing or not.
        Used to update the drag/focus capabilities
        """
        self.canvas.play_overlay.hide_pause(is_playing)
        if CAN_DRAG in self._orig_abilities and self._microscope_view.has_stage():
            # disable/enable move
            if is_playing:
                self.canvas.abilities.add(CAN_DRAG)
            else:
                self.canvas.abilities.discard(CAN_DRAG)
        # check focus ability too
        self._on_stream_change()

    def _on_stream_change(self, streams=None):
        """
        Called whenever the current (playing) stream changes.
        Used to update the focus capability based on the stream
        """
        if CAN_FOCUS not in self._orig_abilities:
            return
        # find out the current playing stream in the view
        for s in self._microscope_view.getStreams():
            if s.should_update.value:
                can_focus = s.focuser is not None
                logging.debug("current stream can focus: %s", can_focus)
                break
        else:
            logging.debug("Found no playing stream")
            can_focus = False

        if can_focus:
            self.canvas.abilities.add(CAN_FOCUS)
        else:
            self.canvas.abilities.discard(CAN_FOCUS)


# TODO: remove once SparcAcquiCanvas is just a normal canvas
class SparcLiveViewport(LiveViewport):
    """
    LiveViewport which support additional modes for the SPARC
    """

    canvas_class = miccanvas.SparcAcquiCanvas


class ARLiveViewport(LiveViewport):
    """
    LiveViewport dedicated to show AR images.
    Never allow to move/zoom, and do not show pause icon if no stream.
    """

    canvas_class = miccanvas.SparcARCanvas

    def __init__(self, *args, **kwargs):
        super(ARLiveViewport, self).__init__(*args, **kwargs)
        # TODO: should be done on the fly by _checkMergeSliderDisplay()
        # change SEM icon to Goal
        if self.bottom_legend:
            self.bottom_legend.bmp_slider_right.SetBitmap(getico_blending_goalBitmap())

    def setView(self, microscope_view, tab_data):
        super(ARLiveViewport, self).setView(microscope_view, tab_data)
        microscope_view.lastUpdate.subscribe(self._on_stream_update, init=True)

    def _on_stream_update(self, _):
        """
        Hide the play icon overlay if no stream are present
        """
        show = len(self._microscope_view.stream_tree) > 0
        self.canvas.play_overlay.show = show


class ARAcquiViewport(ARLiveViewport):
    """
    Same as ARLiveViewport, but without legend
    """
    bottom_legend_class = None


class AngularResolvedViewport(ViewPort):
    """
    Viewport to show the (static) AR images with polar projection
    """

    # Default class
    canvas_class = miccanvas.AngularResolvedCanvas
    bottom_legend_class = None

    def setView(self, microscope_view, tab_data):
        assert(self._microscope_view is None)

        self._microscope_view = microscope_view
        self._tab_data_model = tab_data

        # canvas handles also directly some of the view properties
        self.canvas.setView(microscope_view, tab_data)


class PlotViewport(ViewPort):
    """ Class for displaying plotted data """
    __metaclass__ = ABCMeta

    # Default class
    canvas_class = miccanvas.BarPlotCanvas
    bottom_legend_class = AxisLegend
    left_legend_class = AxisLegend

    def __init__(self, *args, **kwargs):
        super(PlotViewport, self).__init__(*args, **kwargs)
        # We need a local reference to the stream, because if we rely
        # on the reference within the MicroscopeView, it might be replaced
        # before we get an explicit chance to unsubscribe event handlers
        self.stream = None

    def setView(self, view, tab_data):
        """
        Set the microscope view that this viewport is displaying/representing
        *Important*: Should be called only once, at initialisation.

        :param view:(model.View)
        :param tab_data: (model.MicroscopyGUIData)
        """
        # TODO: rename `microscope_view`, since this parameter is a regular view

        # This is a kind of a kludge, as it'd be best to have the viewport
        # created after the microscope view, but they are created independently
        # via XRC.
        assert(self._microscope_view is None)

        # import traceback
        # traceback.print_stack()

        self._microscope_view = view
        self._tab_data_model = tab_data

        # canvas handles also directly some of the view properties
        self.canvas.setView(view, tab_data)

        # Keep an eye on the stream tree, so we can (re)connect when it changes
        # microscope_view.stream_tree.should_update.subscribe(self.connect_stream)
        # FIXME: it shouldn't listen to should_update, but to modifications of
        # the stream tree itself... it just there is nothing to do that.
        view.lastUpdate.subscribe(self.connect_stream)

        view.stream_tree.should_update.subscribe(self._on_stream_play, init=True)
        view.lastUpdate.subscribe(self._on_stream_update, init=True)

    def _on_stream_update(self, _):
        """
        Hide the play icon overlay if no stream are present (or they are all static)
        """
        ss = self._microscope_view.getStreams()
        if len(ss) > 0:
            # Any stream not static?
            show = any(not isinstance(s, StaticStream) for s in ss)
        else:
            show = False
        self.canvas.play_overlay.show = show

    def _on_stream_play(self, is_playing):
        """
        Update the status of the play/pause icon overlay
        """
        self.canvas.play_overlay.hide_pause(is_playing)

    def Refresh(self, *args, **kwargs):
        """
        Refresh the ViewPort while making sure the legends get redrawn as well
        Can be called safely from other threads.
        """
        self.left_legend.Refresh()
        self.bottom_legend.Refresh()
        self.canvas.Refresh()

    def connect_stream(self, _):
        """
        Find the most appropriate stream in the view to be displayed, and make
        sure the display is updated when the stream is updated.
        """
        ss = self._microscope_view.getStreams()
        # Most of the time, there is only one stream, but in some cases, there
        # might be more.
        # TODO: filter based on the type of stream?
        # ss = self.microscope_view.stream_tree.get_streams_by_type(MonochromatorSettingsStream)
        if not ss:
            stream = None
        elif len(ss) > 1:
            # => pick the first one playing
            for s in ss:
                if s.should_update.value:
                    stream = s
                    break
            else:  # no stream playing
                logging.warning("Found %d streams, will pick one randomly", len(ss))
                if self.stream in ss:
                    stream = self.stream  # don't change
                else:
                    stream = ss[0]
        else:
            stream = ss[0]

        if self.stream is stream:
            # logging.debug("not reconnecting to stream as it's already connected")
            return

        # Disconnect the old stream
        if self.stream:
            logging.debug("Disconnecting %s from plotviewport", stream)
            if hasattr(self.stream, 'selected_pixel'):
                self.stream.selected_pixel.unsubscribe(self._on_pixel_select)
            elif hasattr(self.stream, 'image'):
                self.stream.image.unsubscribe(self._on_new_data)

        # Connect the new one
        self.stream = stream
        if stream:
            logging.debug("Connecting %s to plotviewport", stream)
            # Hack: StaticSpectrumStream contain a 2D spectrum in .image, and
            # to get the point spectrum we need to use get_pixel_spectrum() and
            # listen to selected_pixel VA.
            if hasattr(self.stream, 'selected_pixel'):
                self.stream.selected_pixel.subscribe(self._on_pixel_select, init=True)
            elif hasattr(self.stream, 'image'):
                self.stream.image.subscribe(self._on_new_data, init=True)
        else:
            logging.info("No stream to plot found")
            self.clear()  # Remove legend ticks and clear plot

    @abstractmethod
    def _on_new_data(self, data):
        pass

    def _on_pixel_select(self, pixel):
        raise NotImplementedError("This viewport doesn't support streams with .selected_pixel")


class PointSpectrumViewport(PlotViewport):
    """
    Shows the spectrum of a point -> bar plot + legend
    Legend axes are wavelength/intensity.
    """

    def setView(self, view, tab_data):
        super(PointSpectrumViewport, self).setView(view, tab_data)
        wx.CallAfter(self.bottom_legend.SetToolTipString, "Wavelength")
        wx.CallAfter(self.left_legend.SetToolTipString, "Intensity")

    def _on_new_data(self, data):
        """
        Called when a new data is available.
        data (1D DataArray)
        """
        if data.size:
            # TODO: factorize with get_spectrum_range() for static stream?
            try:
                spectrum_range = spectrum.get_wavelength_per_pixel(data)
                unit_x = "m"
            except (ValueError, KeyError):
                # useless polynomial => just show pixels values (ex: -50 -> +50 px)
                max_bw = data.shape[0] // 2
                min_bw = (max_bw - data.shape[0]) + 1
                spectrum_range = range(min_bw, max_bw + 1)
                unit_x = "px"

            self.canvas.set_1d_data(spectrum_range, data, unit_x)

            self.bottom_legend.unit = unit_x
            self.bottom_legend.range = (spectrum_range[0], spectrum_range[-1])
            self.left_legend.range = (min(data), max(data))
            # For testing
            # self.left_legend.range = (min(data) + random.randint(0, 100), max(data) + random.randint(-100, 100))
        else:
            self.clear()
        self.Refresh()

    def _on_pixel_select(self, pixel):
        """ Pixel selection event handler """
        if pixel == (None, None):
            # TODO: handle more graciously when pixel is unselected?
            logging.debug("No pixel selected")
            # Remove legend ticks and clear plot
            self.clear()
            return
        elif self.stream is None:
            logging.warning("No Spectrum Stream present!")
            return

        data = self.stream.get_pixel_spectrum()
        spectrum_range = self.stream.get_spectrum_range()
        unit_x = self.stream.spectrumBandwidth.unit

        self.canvas.set_1d_data(spectrum_range, data, unit_x)

        self.bottom_legend.unit = unit_x
        self.bottom_legend.range = (spectrum_range[0], spectrum_range[-1])
        self.left_legend.range = (min(data), max(data))

        self.Refresh()


class ChronographViewport(PlotViewport):
    """
    Shows the chronograph of a 0D detector reading -> bar plot + legend
    Legend axes are time/intensity.
    """

    def __init__(self, *args, **kwargs):
        super(ChronographViewport, self).__init__(*args, **kwargs)
        self.canvas.markline_overlay.hide_x_label()

    def setView(self, view, tab_data):
        super(ChronographViewport, self).setView(view, tab_data)
        wx.CallAfter(self.bottom_legend.SetToolTipString, "Time (s)")
        wx.CallAfter(self.left_legend.SetToolTipString, "Count per second")

    def _on_new_data(self, data):
        if data.size:
            unit_x = 's'

            x = data.metadata[model.MD_ACQ_DATE]
            y = data
            range_x = (min(x[0], -self.stream.windowPeriod.value), x[-1])
            range_y = (0, float(max(data)))  # float() to avoid numpy arrays

            self.canvas.set_data(zip(x, y), unit_x, range_x=range_x, range_y=range_y)

            self.bottom_legend.unit = unit_x
            self.bottom_legend.range = range_x
            self.left_legend.range = range_y

        else:
            self.clear()
        self.Refresh()


class SpatialSpectrumViewport(ViewPort):
    """
    A viewport for showing 1D spectum: an image with wavelength horizontally and
    space vertically.
    """
    # FIXME: This class shares a lot with PlotViewport, see what can be merged

    canvas_class = miccanvas.TwoDPlotCanvas
    bottom_legend_class = AxisLegend
    left_legend_class = AxisLegend

    def __init__(self, *args, **kwargs):
        """Note: The MicroscopeViewport is not fully initialised until setView()
        has been called.
        """
        # Call parent constructor at the end, because it needs the legend panel
        super(SpatialSpectrumViewport, self).__init__(*args, **kwargs)
        self.stream = None
        self.current_line = None

        self.canvas.markline_overlay.v_pos.subscribe(self.on_spectrum_motion)

    def on_spectrum_motion(self, vpos):

        if vpos:
            self.canvas.markline_overlay.x_label = units.readable_str(
                self.bottom_legend.pixel_to_value(vpos[0]),
                self.bottom_legend.unit,
                3
            )
            self.canvas.markline_overlay.y_label = units.readable_str(
                self.left_legend.pixel_to_value(vpos[1]),
                self.left_legend.unit,
                3
            )
            rat = self.left_legend.pixel_to_ratio(vpos[1])
            line_pixels = rasterize_line(*self.current_line)
            self.stream.selected_pixel.value = line_pixels[int(len(line_pixels) * rat)]

    def Refresh(self, *args, **kwargs):
        """
        Refresh the ViewPort while making sure the legends get redrawn as well
        Can be called safely from other threads
        """
        self.left_legend.Refresh()
        self.bottom_legend.Refresh()
        # Note: this is not thread safe, so would need to be in a CallAfter()
        # super(SpatialSpectrumViewport, self).Refresh(*args, **kwargs)
        wx.CallAfter(self.canvas.update_drawing)

    def setView(self, microscope_view, tab_data):
        """
        Set the microscope view that this viewport is displaying/representing
        *Important*: Should be called only once, at initialisation.

        :param microscope_view:(model.View)
        :param tab_data: (model.MicroscopyGUIData)
        """
        # TODO: rename `microscope_view`, since this parameter is a regular view

        # This is a kind of a kludge, as it'd be best to have the viewport
        # created after the microscope view, but they are created independently
        # via XRC.
        assert(self._microscope_view is None)

        # import traceback
        # traceback.print_stack()

        self._microscope_view = microscope_view
        self._tab_data_model = tab_data

        # canvas handles also directly some of the view properties
        self.canvas.setView(microscope_view, tab_data)

        # Keep an eye on the stream tree, so we can (re)connect when it changes
        # microscope_view.stream_tree.should_update.subscribe(self.connect_stream)
        # FIXME: it shouldn't listen to should_update, but to modifications of
        # the stream tree itself... it just there is nothing to do that.
        microscope_view.lastUpdate.subscribe(self.connect_stream)

        wx.CallAfter(self.bottom_legend.SetToolTipString, "Wavelength")
        wx.CallAfter(self.left_legend.SetToolTipString, "Distance from origin")

    def connect_stream(self, _=None):
        """ This method will connect this ViewPort to the Spectrum Stream so it
        it can react to spectrum pixel selection.
        """
        ss = self.microscope_view.stream_tree.spectrum_streams
        if self.stream in ss:
            logging.debug("not reconnecting to stream as it's already connected")
            return

        # There should be exactly one Spectrum stream. In the future there
        # might be scenarios where there are more than one.
        if not ss:
            self.stream = None
            logging.info("No spectrum streams found")
            self.clear()  # Remove legend ticks and clear image
            return
        elif len(ss) > 1:
            logging.warning("Found %d spectrum streams, will pick one randomly", len(ss))

        self.stream = ss[0]
        self.stream.selected_line.subscribe(self._on_line_select, init=True)
        self.stream.selected_pixel.subscribe(self._on_pixel_select)

    def _on_pixel_select(self, pixel):
        """ Clear the marking line when the selected pixel is cleared """
        if None in pixel:
            self.canvas.markline_overlay.clear_labels()

    def _on_line_select(self, line):
        """ Line selection event handler """

        if (None, None) in line:
            logging.debug("Line is not (fully) selected")
            self.clear()
            self.current_line = None
            return
        elif self.stream is None:
            logging.warning("No Spectrum Stream present!")
            return

        data = self.stream.get_line_spectrum()
        self.current_line = line

        if data is not None:
            spectrum_range = self.stream.get_spectrum_range()
            line_length = data.shape[0] * data.metadata[model.MD_PIXEL_SIZE][1]

            unit_x = self.stream.spectrumBandwidth.unit
            self.bottom_legend.unit = unit_x
            self.bottom_legend.range = (spectrum_range[0], spectrum_range[-1])
            self.left_legend.unit = 'm'
            self.left_legend.range = (0, line_length)

            self.canvas.set_2d_data(data)
        else:
            logging.warn("No data to display for the selected line!")

        self.Refresh()
