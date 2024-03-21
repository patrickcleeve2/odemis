
# -*- coding: utf-8 -*-

"""
@author Patrick Cleeve

Copyright © 2023, Delmic

Handles the controls for correlating two (or more) streams together.

This file is part of Odemis.

Odemis is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License version 2 as published by the Free
Software Foundation.

Odemis is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
Odemis. If not, see http://www.gnu.org/licenses/.

"""
import logging
import math
import wx
# IMPORTANT: wx.html needs to be imported for the HTMLWindow defined in the XRC
# file to be correctly identified. See: http://trac.wxwidgets.org/ticket/3626
# This is not related to any particular wxPython version and is most likely permanent.
import wx.html
import odemis.acq.stream as acqstream
import odemis.gui.model as guimod
from odemis import model
from odemis.acq.stream import StaticStream, StaticFluoStream, StaticSEMStream
from odemis.gui.util import call_in_wx_main
from odemis.gui.cont.tabs.localization_tab import LocalizationTab

# TODO: move to more approprate location
def update_image_in_views(s: StaticStream, views: list) -> None:
    """Force update the static stream in the selected views
    :param s: (StaticStream) the static stream to update
    :param views: (list[View]) the list of views to update"""
    v: guimod.View
    for v in views:
        for sp in v.stream_tree.getProjections():  # stream or projection
            if isinstance(sp, acqstream.DataProjection):
                st = sp.stream
            else:
                st = sp
            if st is s:
                sp._shouldUpdateImage()


class CorrelationController(object):

    def __init__(self, tab_data, panel, tab, viewports):
        """
        :param tab_data: (MicroscopyGUIData) the representation of the microscope GUI
        :param panel: (wx._windows.Panel) the panel containing the UI controls
        :param tab: (Tab) the tab which should show the data
        :param viewports: (ViewPorts) the tab view ports
        """

        self._tab_data_model = tab_data
        self._main_data_model = tab_data.main
        self._panel = panel
        self._tab = tab
        self._viewports = viewports

        # disable if no streams are present
        self._tab_data_model.streams.subscribe(self._on_correlation_streams_change, init=True)
        self._tab_data_model.streams.subscribe(self._update_correlation_cmb, init=True)
        self._tab_data_model.streams.subscribe(self.add_to_localization_tab, init=False)
        self._tab_data_model.streams.subscribe(self._update_point_correlation_cmb, init=True)

        # connect the correlation streams to the tab data
        self._panel.cmb_correlation_stream.Bind(wx.EVT_COMBOBOX, self._on_selected_stream_change)

        # reference frames
        self.ref_frame_names = ["METEOR", "FIBSEM", "No Reference"] # for display purposes
        self.ref_stream_map = {0: StaticFluoStream,  # METEOR
                               1: StaticSEMStream,  # FIBSEM
                               2: type(None)} # match everything except None (No reference type)
        self.ref_stream = self.ref_stream_map[0]

        # connect the reference stream to the tab data
        self._panel.cmb_correlation_reference.Append(self.ref_frame_names)
        self._panel.cmb_correlation_reference.SetSelection(0)
        self._panel.cmb_correlation_reference.Bind(wx.EVT_COMBOBOX, self._on_reference_stream_change)

        # reset correlation data
        self._panel.btn_reset_correlation.Bind(wx.EVT_BUTTON, self.reset_correlation_pressed)

        # enable correlation controls
        self._panel.ctrl_enable_correlation.SetValue(True) # enable by default
        self._panel.ctrl_auto_resize_view.SetValue(False) # disable auto resize by default
        self._panel.ctrl_enable_correlation.SetToolTip("Enable/Disable correlation controls")
        self._panel.ctrl_auto_resize_view.SetToolTip("Automatically resize correlation overlay viewports after moving the streams.")

        # point correlation
        self.sem_features = model.ListVA()
        self.flm_features = model.ListVA()

        self.sem_features.subscribe(self._on_features_changes, init=True)
        self.flm_features.subscribe(self._on_features_changes, init=True)

        # self.points_overlay = None
        # self.add_world_overlay(self.points_overlay)
        # self.points_overlay.active.value = True

        # bind mouse and keyboard events for correlation controls
        for vp in self._viewports:
            vp.canvas.Bind(wx.EVT_CHAR, self.on_char)
            vp.canvas.Bind(wx.EVT_LEFT_DOWN, self.on_mouse_down)

        # localization tab
        self.localization_tab: LocalizationTab = None

    @call_in_wx_main
    def _update_correlation_cmb(self, streams: list) -> None:
        """update the correlation combo box with the available streams
        :param streams: (list[StaticStream]) the streams to add to the combo box"""
        # keep the combobox in sync with streams
        self._panel.cmb_correlation_stream.Clear()
        for s in streams:
            if isinstance(s, self.ref_stream):
                continue # cant move fluo streams
            self._panel.cmb_correlation_stream.Append(s.name.value, s)

        # select the first stream, if available and nothing is selected
        if (self._panel.cmb_correlation_stream.GetCount() > 0
            and self._panel.cmb_correlation_stream.GetSelection() == wx.NOT_FOUND):
            self._panel.cmb_correlation_stream.SetSelection(0) # this doesn't trigger selection event

            # trigger the event manually, to automatically select the first stream
            if self._tab_data_model.selected_stream.value is None:
                logging.debug(f"Forcing selected stream change to first stream")
                self._on_selected_stream_change(None)

    @call_in_wx_main
    def _update_point_correlation_cmb(self, streams: list) -> None:
        """update the point correlation comboboxes with the available streams
        :param streams (list[StaticStream]) the streams to add to the combo box"""

        self._panel.cmb_point_correlation_reference.Clear()
        self._panel.cmb_point_correlation_moving.Clear()
        for s in streams:
            self._panel.cmb_point_correlation_reference.Append(s.name.value, s)
            self._panel.cmb_point_correlation_moving.Append(s.name.value, s)
        
        if self._panel.cmb_point_correlation_reference.GetCount() > 0:
            self._panel.cmb_point_correlation_reference.SetSelection(0)

        if self._panel.cmb_point_correlation_moving.GetCount() > 0:
            self._panel.cmb_point_correlation_moving.SetSelection(0)    

    @call_in_wx_main
    def _on_correlation_streams_change(self, streams: list) -> None:
        """hide/show the correlation panel if there are no streams
        :param streams: (list[StaticStream]) the streams in the correlation tab"""
        visible = len(streams) != 0
        self._panel.fp_meteor_correlation.Show(visible)
        # reset selected stream
        if not visible:
            self._panel.cmb_correlation_stream.Clear()
            self._tab_data_model.selected_stream.value = None

    @call_in_wx_main
    def _on_selected_stream_change(self, evt: wx.Event) -> None:
        """change the selected stream to the one selected in the combo box
        :param evt: (wx.Event) the event"""
        idx = self._panel.cmb_correlation_stream.GetSelection()
        self._tab_data_model.selected_stream.value = self._panel.cmb_correlation_stream.GetClientData(idx)
        logging.debug(f"Selected Stream Changed to {idx}: {self._tab_data_model.selected_stream.value.name.value}")

    @call_in_wx_main
    def _on_reference_stream_change(self, evt: wx.Event) -> None:
        """change the reference stream to the one selected in the combo box
        :param evt: (wx.Event) the event"""
        idx = self._panel.cmb_correlation_reference.GetSelection()
        self.ref_stream = self.ref_stream_map[idx]
        logging.debug(f"Reference Frame: {self.ref_frame_names[idx]} - Fixed stream: {idx}, {self.ref_stream}")

        # refresh combobox
        self._update_correlation_cmb(self._tab_data_model.streams.value)

    def reset_correlation_pressed(self, evt: wx.Event) -> None:
        """"Reset the correlation data for the selected stream, and re-draw
        :param evt: (wx.Event) the event"""
        s = self._tab_data_model.selected_stream.value
        self._reset_stream_correlation_data(s)
        update_image_in_views(s, self._tab_data_model.views.value)
        self.fit_correlation_views_to_content(force=True)

    def _reset_stream_correlation_data(self, s: StaticStream) -> None:
        """reset the stream position to the original position / rotation / scale
        :param s: (StaticStream) the stream to reset"""
        if s is not None and s.raw:
            s.raw[0].metadata[model.MD_POS_COR] = (0, 0)
            s.raw[0].metadata[model.MD_ROTATION_COR] = 0
            s.raw[0].metadata[model.MD_PIXEL_SIZE_COR] = (1, 1)

    def add_streams(self, streams: list) -> None:
        """add streams to the correlation tab
        :param streams: (list[StaticStream]) the streams to add"""

        # NOTE: we only add streams if they are not already in the correlation tab
        # we will not remove streams from the correlation tab from outside it,
        # as the user may still want to correlate them,
        # even if they have been removed from another tab, e.g. 'localization'

        # add streams to correlation tab
        logging.debug(f"Adding {len(streams)} streams to correlation tab {streams}")
        for s in streams:

            # skip existing streams, live streams
            if s in self._tab_data_model.streams.value or not isinstance(s, StaticStream):
                continue

            # reset the stream correlation data, and add to correlation streams
            self._reset_stream_correlation_data(s)
            self._tab_data_model.streams.value.append(s)

            # add stream to streambar
            sc = self._tab.streambar_controller.addStream(s, add_to_view=True, play=False)
            sc.stream_panel.show_remove_btn(True)

    def add_to_localization_tab(self, streams: list) -> None:
        """add streams to the localization tab
        :param streams: (list[StaticStream]) the streams to add"""
        # NOTE: we only add streams if they are not already in the localization tab
        # we will not remove streams from the localization tab from outside it,
        # as the user may still want to localize them,
        # even if they have been removed from another tab, e.g. 'correlation'
        if self.localization_tab is None:
            self.localization_tab: LocalizationTab  = self._main_data_model.getTabByName("cryosecom-localization")

        logging.debug(f"Adding {len(streams)} streams to localization tab {streams}")
        for s in streams:
            # add stream to localizations tab
            if s not in self.localization_tab.tab_data_model.streams.value:
                self.localization_tab.tab_data_model.overviewStreams.value.append(s)
                self.localization_tab.tab_data_model.streams.value.insert(0, s)
                self.localization_tab._acquired_stream_controller.showOverviewStream(s)


    def clear_streams(self) -> None:
        """clears streams from the correlation tab"""
        self._tab.streambar_controller.clear()

    def correlation_enabled(self) -> bool:
        """return if correlation controls are enabled and a stream is selected"""
        logging.debug(f"correlation enabled: {self._panel.ctrl_enable_correlation.IsChecked()}")
        logging.debug(f"selected stream: {self._tab_data_model.selected_stream.value}")
        return (self._panel.ctrl_enable_correlation.IsChecked() and
                self._tab_data_model.selected_stream.value is not None)

    def on_mouse_down(self, evt: wx.Event) -> None:
        """handle mouse down events
        :param evt: (wx.Event) the event"""
        active_canvas = evt.GetEventObject()
        logging.debug(f"mouse down event, canvas: {active_canvas}")

        # check if shift is pressed, and if a stream is selected
        if evt.ShiftDown() and self.correlation_enabled():

            # get the position of the mouse, convert to physical position
            pos = evt.GetPosition()
            p_pos = active_canvas.view_to_phys(pos, active_canvas.get_half_buffer_size())
            logging.debug(f"shift pressed, mouse_pos: {pos}, phys_pos: {p_pos}")

            # move selected stream to position
            self._move_stream_to_pos(p_pos)
        
        elif evt.ControlDown():

            # add a feature to canvas
            pos = evt.GetPosition()

            p_pos = active_canvas.view_to_phys(pos, active_canvas.get_half_buffer_size())

            logging.info(f"FEATURE POSITION: {p_pos}")
            
            # add a feature to the canvas
            if active_canvas == self._panel.vp_correlation_tl.canvas:
                self._add_feature(p_pos, "FLM")
            elif active_canvas == self._panel.vp_correlation_tr.canvas:
                self._add_feature(p_pos, "SEM")
            else:
                logging.info("No feature added, wrong canvas")
        else:
            logging.debug("invalid correlation event, passing event to canvas")
            active_canvas.on_left_down(evt)       # super event passthrough

    def on_char(self, evt: wx.Event) -> None:
        """handle key presses
        :param evt: (wx.Event) the event"""

        # event data
        key = evt.GetKeyCode()
        shift_mod = evt.ShiftDown()

        # pass through event, if not a valid correlation key or enabled
        valid_keys = [wx.WXK_LEFT, wx.WXK_RIGHT, wx.WXK_UP, wx.WXK_DOWN, wx.WXK_SPACE, wx.WXK_DELETE]
        if key not in valid_keys or not self.correlation_enabled():
            evt.Skip()
            return

        ### CONTROLS ##############################
        # SHIFT + LEFT CLICK -> MOVE_TO_POSITION
        # LEFT, RIGHT -> TRANSLATION X
        # UP, DOWN -> TRANSLATION Y
        # SHIFT + LEFT, RIGHT -> ROTATION
        # SHIFT + UP, DOWN -> SCALE
        ###########################################
        dx, dy, dr, dpx  = 0, 0, 0, 0

        # correlation control modifiers
        if shift_mod:
            dr = math.radians(self._panel.dr_step_cntrl.GetValue())
            dpx = self._panel.dpx_step_cntrl.GetValue() / 100
        else:
            dx = self._panel.dx_step_cntrl.GetValue()
            dy = self._panel.dy_step_cntrl.GetValue()

        logging.debug(f"key: {key}, shift: {shift_mod}")
        logging.debug(f"dx: {dx}, dy: {dy}, dr: {dr}, dpx: {dpx}")

        if key == wx.WXK_LEFT:
            self._move_stream(-dx, 0, -dr, 0)
        elif key == wx.WXK_RIGHT:
            self._move_stream(dx, 0, dr, 0)
        elif key == wx.WXK_UP:
            self._move_stream(0, dy, 0, dpx)
        elif key == wx.WXK_DOWN:
            self._move_stream(0, -dy, 0, -dpx)


        
        if key == wx.WXK_SPACE:
            # test dialog
            # if self.sem_stream is None:       
            #     box = wx.MessageDialog(self.main_frame,
            #             "No stream is present, so it's not possible to modify the alignment.",
            #             "No stream", wx.OK | wx.ICON_STOP)
            #     box.ShowModal()
            #     box.Destroy()
            #     return
            from odemis.acq.feature import get_features_dict
            logging.info(f"SEM Features: {get_features_dict(self.sem_features.value)}")
            logging.info(f"FLM Features: {get_features_dict(self.flm_features.value)}")


        if key == wx.WXK_DELETE:
            active_canvas = evt.GetEventObject()
            # delete all features from canvas
            if active_canvas == self._panel.vp_correlation_tl.canvas:
                feature_type = "SEM"
            elif active_canvas == self._panel.vp_correlation_tr.canvas:
                feature_type = "FLM"
            else:
                return
            self._delete_features(feature_type)

    def _move_stream(self, dx: float, dy: float , dr: float = 0, dpx: float = 0) -> None:
        """move the selected stream by the specified amount. the stream is forced
        to update in the views.
        :param dx: (float) the change in x translation (metres)
        :param dy: (float) the change y translation (metres)
        :param dr: (float) the change in rotation (radians)
        :param dpx: (float) the change in scale of the pixelsize (percentage)
        """
        if self._tab_data_model.selected_stream.value is None:
            return

        s = self._tab_data_model.selected_stream.value

        logging.debug(f"move stream {s.name.value}: {dx}, {dy}, {dr}, {dpx}")

        # translation
        p = s.raw[0].metadata[model.MD_POS_COR]
        if len(p) == 2:
            x, y = p
            s.raw[0].metadata[model.MD_POS_COR] = (x - dx, y - dy) # correlation direction is reversed
        else:
            x, y, z = p
            s.raw[0].metadata[model.MD_POS_COR] = (x - dx, y - dy, z)

        # rotation
        rotation = s.raw[0].metadata.get(model.MD_ROTATION_COR, 0)
        s.raw[0].metadata[model.MD_ROTATION_COR] = (rotation + dr)

        # scale (pixel size)
        scalecor = s.raw[0].metadata.get(model.MD_PIXEL_SIZE_COR, (1, 1))
        s.raw[0].metadata[model.MD_PIXEL_SIZE_COR] = (scalecor[0] + dpx, scalecor[1] + dpx)

        # TODO: split x, y scale, add shear?

        # update the image in the views
        update_image_in_views(s, self._tab_data_model.views.value)

        # fit viewports to content, as they have moved
        self.fit_correlation_views_to_content()

        # update the localization tab
        self.update_localization_tab(s)

    def fit_correlation_views_to_content(self, force: bool = False) -> None:
        # TODO: be more selective about which viewports to fit
        # fit the source viewports to the content, as the image may have moved
        self._panel.vp_correlation_tl.canvas.fit_view_to_content()
        self._panel.vp_correlation_tr.canvas.fit_view_to_content()

        # fit the overlay viewports to the content (optional)
        auto_resize = self._panel.ctrl_auto_resize_view.IsChecked()
        if auto_resize or force:
            self._panel.vp_correlation_bl.canvas.fit_view_to_content()
            self._panel.vp_correlation_br.canvas.fit_view_to_content()

    def update_localization_tab(self, s: StaticStream) -> None:
        # TODO: change this to callback?
        # also update the localization tab
        if s in self.localization_tab.tab_data_model.streams.value:

            logging.debug(f"Updating localisation stream: {s.name.value}")

            # get stream in localization tab
            idx = self.localization_tab.tab_data_model.streams.value.index(s)
            sl = self.localization_tab.tab_data_model.streams.value[idx]

            # match cor metadata
            sl.raw[0].metadata[model.MD_POS_COR] = s.raw[0].metadata[model.MD_POS_COR]
            sl.raw[0].metadata[model.MD_ROTATION_COR] = s.raw[0].metadata[model.MD_ROTATION_COR]
            sl.raw[0].metadata[model.MD_PIXEL_SIZE_COR] = s.raw[0].metadata[model.MD_PIXEL_SIZE_COR]

            update_image_in_views(s, self.localization_tab.tab_data_model.views.value)

    def _move_stream_to_pos(self, pos: tuple) -> None:
        """move the selected stream to the position pos
        :param pos: (tuple) the realspace position to move the stream to (metres)"""
        # the difference between the clicked position, and the position in metadata
        # is the offset (be careful because correlation is sign flipped)
        s = self._tab_data_model.selected_stream.value

        # the cur_pos is the realspace position of the image
        p = s.raw[0].metadata[model.MD_POS]
        x, y = p[:2]  # x, y positions only (ignore z)

        # the correlation pos is the change in position in the viewer
        cx, cy = s.raw[0].metadata[model.MD_POS_COR]

        # the new position is the difference between the clicked position and the realspace position
        # i.e. the correlation position. However, already have an existing correlation position, so we need to
        # find how much we need to modify this by.
        nx = -(pos[0] - x)
        ny = -(pos[1] - y)

        # the offset is the difference between the current cor position and the new position
        # that is how much we need to move the current correlation position by.
        dx = cx - nx
        dy = cy - ny

        logging.debug(f"cur_pos: {x}, {y}, cor_pos: {cx}, {cy}, new_pos: {nx}, {ny}, offset_pos: {dx}, {dy}")

        # move the stream using the correlation position offset
        self._move_stream(dx=dx, dy=dy, dr=0, dpx=0)


    def _add_feature(self, pos: tuple, feature_type: str):
        
        from odemis.acq.feature import CryoFeature
        
        _num = len(self.sem_features.value) if feature_type == "SEM" else len(self.flm_features.value)

        f = CryoFeature(name=f"{feature_type} Feature {_num}", 
                            x=pos[0], y=pos[1], z=0)
        if feature_type == "SEM":
            self.sem_features.value.append(f)
        elif feature_type == "FLM":
            self.flm_features.value.append(f)
        else:
            logging.info("No feature added, wrong feature type")
            return
        
        self._tab_data_model.main.features.value.append(f)

    def _delete_features(self, feature_type: str):

        # remove all features from canvas
        if feature_type == "SEM":
            while len(self.sem_features.value) > 0:
                self.sem_features.value.pop()
        elif feature_type == "FLM":
            while len(self.flm_features.value) > 0:
                self.flm_features.value.pop()
        else:
            logging.info("No feature deleted, wrong feature type")
        
        self._tab_data_model.main.features.value = []
    
    def _on_features_changes(self, features):
        from odemis.acq.feature import get_features_dict

        # check if features is empty
        if not features:
            return
        
        logging.info(f"SEM Features: {get_features_dict(self.sem_features.value)}")
        logging.info(f"FLM Features: {get_features_dict(self.flm_features.value)}")

        logging.info(f"SEM FEATURES: {len(self.sem_features.value)}")
        logging.info(f"FLM FEATURES: {len(self.flm_features.value)}")

        MIN_REQ_FEATURES = 4
        if len(self.sem_features.value) >= MIN_REQ_FEATURES and len(self.flm_features.value) >= MIN_REQ_FEATURES:
        # and equal length
            if len(self.sem_features.value) == len(self.flm_features.value):
                logging.info(f"FEATURES EQUAL LENGTH: DO TRANSFORMATION")
                
                import numpy as np
                # TODO: implement

                # calculate the order of points

                sem_points = [f.pos.value for f in self.sem_features.value]
                flm_points = [f.pos.value for f in self.flm_features.value]


                # assume always flm -> sem
                # calculate transformation matrix
                src = np.array(flm_points)
                dst = np.array(sem_points)

                # ODEMIS TRANSFORM
                from odemis.util.transform import AffineTransform

                # convert src, dst to 2d points
                src = [f[:2] for f in src]
                dst = [f[:2] for f in dst]

                affine1 = AffineTransform.from_pointset(src, dst)
                affine2 = AffineTransform.from_pointset(dst, src)
                
                logging.info(f"-----------------------")
                logging.info(f"MATRIX 1: {affine1.matrix}")
                logging.info(f"TRANSLATION 1: {affine1.translation}")
                logging.info(f"-----------------------")
                logging.info(f"MATRIX 2: {affine2.matrix}")
                logging.info(f"TRANSLATION 2: {affine2.translation}")
                logging.info(f"-----------------------")
                
                # apply transformation to the selected stream



                # TODO: to matrix multiplication
                from odemis.util.conversion import get_img_transformation_md, get_img_transformation_matrix
                from odemis.util.transform import alt_transformation_matrix_to_implicit


                flm_stream = None
                sem_stream = None

                # get the position of the flm stream
                for stream in self._tab_data_model.streams.value:
                    if isinstance(stream, StaticFluoStream):
                        flm_stream = stream
                        break


                # get the position of the sem stream
                for stream in self._tab_data_model.streams.value:
                    if isinstance(stream, StaticSEMStream):
                        sem_stream = stream
                        break
                

                # TODO: implement user selection of reference, moving stream
                # idx = self._panel.cmb_point_correlation_reference.GetSelection()
                # self._tab_data_model.selected_stream.value = self._panel.cmb_point_correlation_reference.GetClientData(idx)
                # logging.debug(f"Selected Stream Changed to {idx}: {self._tab_data_model.selected_stream.value.name.value}")


                if flm_stream is None or sem_stream is None:
                    logging.info(f"NO STREAMS")
                    return

                # assume always flm -> sem (for now)
                #         
                #REF https://github.com/jni/affinder/blob/main/src/affinder/affinder.py#L74

                import copy


                def get_tsr(mat):
                    # Translation is the last column
                    tx = mat[0, 2]
                    ty = mat[1, 2]
                    translation = np.array([tx, ty])

                    # Scale is the norm of the first two columns
                    sx = np.linalg.norm(mat[:2, 0])
                    sy = np.linalg.norm(mat[:2, 1])
                    scale2 = np.array([sx, sy])

                    # Compute rotation
                    rotation = np.arctan2(mat[1, 0]/sx, mat[0, 0]/sx)

                    logging.info(f"Translation {translation}")
                    logging.info(f"Scale {scale2}")
                    logging.info(f"Rotation {rotation}")

                    return translation, scale2, rotation
                
                USE_REF = False

                # V1
                if USE_REF:
    
                    ref_trans = sem_stream.raw[0].metadata[model.MD_POS]

                    ref_mat = np.eye(3) # reference is identity. 
                    # ref_mat[:2, :2] = ref_m
                    # ref_mat[:2, 2] = ref_trans

                    logging.info(f"REF MAT: {ref_mat}")

                    # Create a 3x3 identity matrix
                    mat = np.eye(3)
                    mat[:2, :2] = affine2.matrix
                    mat[:2, 2] = affine2.translation
                    
                    # transform from ref to moving
                    mat = ref_mat @ mat

                    translation, scale2, rotation = get_tsr(mat)
                    tx, ty = translation
                    # set MD_PIXEL_SIZE_COR to 1/scale2
                    flm_stream.raw[0].metadata[model.MD_PIXEL_SIZE_COR] = (1/scale2[0], 1/scale2[1])

                    # set MD_POS_COR to translation
                    sem_pos = sem_stream.raw[0].metadata[model.MD_POS]
                    flm_stream.raw[0].metadata[model.MD_POS_COR] = (tx, ty)

                    # set MD_ROTATION_COR to rotation
                    flm_stream.raw[0].metadata[model.MD_ROTATION_COR] = rotation

                    # classically, we would just apply the transformation to the reference frame 
                    # and then apply the inverse to the moving frame.
                    # but odemis uses metadata, not the affine to transform the image
                    # so we need to convert the affine to metadata (t, r, s) and then apply them as corrections
                    # to the moving frame -> MD_POS_COR, MD_ROTATION_COR, MD_PIXEL_SIZE_COR
                    # this should translation, rotate, and scale the moving frame to match the reference frame

                    # alternatively, we can just set the metadata to the transformation matrix directly
                    # MD_POS, MD_ROTATION, MD_PIXEL_SIZE




                ####### V2
                if not USE_REF:
                    mat = np.eye(3)
                    mat[:2, :2] = affine1.matrix
                    mat[:2, 2] = affine1.translation

                    translation, scale2, rotation = get_tsr(mat)

                    # copy initial flm pixel size
                    offset_stream = sem_stream
                    pixel_size  = copy.deepcopy(offset_stream.raw[0].metadata[model.MD_PIXEL_SIZE])
                    new_pixel_size = (pixel_size[0] * scale2[0], pixel_size[1] * scale2[1])
                    transf_md = copy.deepcopy(get_img_transformation_md(mat, flm_stream.raw[0], sem_stream.raw[0]))

                    # subtract half the width / height of the sem image
                    shape = offset_stream.raw[0].shape

                    pos = transf_md[model.MD_POS]

                    # subtract half the width / height of the moving? image
                    pos2 = (pos[0] + (shape[1] / 2 * new_pixel_size[0]), 
                                            pos[1] - (shape[0] / 2 * new_pixel_size[1]))
                
                    if len(pos) == 3:
                        pos2 = (pos2[0], pos2[1], pos[2])

                    logging.info(f"TRANSFORMATION MD: {transf_md}")
                    logging.info(f"TRANSLATION: {translation}")

                    # USING INVERSE TRANSFORMATION
                    # flm_stream.raw[0].metadata[model.MD_POS] = pos2
                    # flm_stream.raw[0].metadata[model.MD_PIXEL_SIZE_COR] = (1/scale2[0], 1/scale2[1])
                    # flm_stream.raw[0].metadata[model.MD_ROTATION_COR] = -transf_md[model.MD_ROTATION]


                    # # calculate the position correction
                    # x_factor, y_factor = np.cos(transf_md[model.MD_ROTATION]), np.sin(transf_md[model.MD_ROTATION])
                    # x_offset = (shape[1] / 2 * new_pixel_size[0]) * x_factor - (shape[0] / 2 * new_pixel_size[1]) * y_factor
                    # y_offset = (shape[1] / 2 * new_pixel_size[0]) * y_factor + (shape[0] / 2 * new_pixel_size[1]) * x_factor

                    # logging.info(f"X OFFSET: {x_offset}, Y OFFSET: {y_offset}")

                    x2, y2 = (shape[1] / 2 * new_pixel_size[0]), (shape[0] / 2 * new_pixel_size[1])
                    logging.info(f"X2: {x2}, Y2: {y2}")

                    # V3
                    flm_stream.raw[0].metadata[model.MD_POS] = pos
                    flm_stream.raw[0].metadata[model.MD_POS_COR] = -x2, y2
                    flm_stream.raw[0].metadata[model.MD_PIXEL_SIZE_COR] = (scale2[0], scale2[1])
                    flm_stream.raw[0].metadata[model.MD_ROTATION_COR] = transf_md[model.MD_ROTATION]


                    # NOTES:
                    # if i just use the pos (translation) directly, the image always ends up in the top left corner of the sem image
                    # it is almost exactly hafl the image width / height away, but not exactly
                    # but not exactly.
                    # im not sure why the transform works like that. Is this a side effect of get_img_transformation_md?
                    # do i just need to add an additional offset to the translation? why? 
                    # scale and rotation are acccurate.
                    # im very confused



                    # add a new flm stream, representing the transformation
                    # flm_transformed = StaticFluoStream(name="flm-transformed", raw=flm_stream.raw[0])
                    # flm_transformed.raw[0].metadata.update(transf_md)
                    # dont match pixel size (scale not working atm)
                    # flm_transformed.raw[0].metadata[model.MD_PIXEL_SIZE] = new_pixel_size

                    # add the new stream to the tab
                    # self.add_streams([flm_transformed])

                    # hide the original flm stream in all views
                    # for v in self._tab_data_model.views.value:
                    #     v.stream_tree.hideStream(flm_stream)
                
                # update the image in the views
                update_image_in_views(flm_stream, self._tab_data_model.views.value)

                # fit the source viewports to the content, as the image may have moved
                self._panel.vp_correlation_tl.canvas.fit_view_to_content()
                self._panel.vp_correlation_tr.canvas.fit_view_to_content()
                self._panel.vp_correlation_bl.canvas.fit_view_to_content()
                self._panel.vp_correlation_br.canvas.fit_view_to_content()
                    
            else:
                logging.info(f"NOT EQUAL LENGTH")
        else:
            logging.info(f"NOT ENOUGH FEATURES")
            