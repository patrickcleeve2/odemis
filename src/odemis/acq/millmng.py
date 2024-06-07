"""
1.Check milling with one preset above the specified lamella
2.Repeat the process to do milling below the lamella

Assumptions:
- lamella is at centre of the FOV

"""
import logging
import math
import os
import threading
import time
from concurrent import futures
from concurrent.futures._base import CANCELLED, FINISHED, RUNNING, CancelledError

import numpy
import yaml

from odemis import model
from odemis.acq import move
from odemis.acq.acqmng import acquire
from odemis.acq.drift import AnchoredEstimator
from odemis.acq.feature import CryoFeature
from odemis.acq.move import MILLING, MicroscopePostureManager
from odemis.acq.orsay_milling import mill_rectangle
from odemis.acq.stitching._tiledacq import MOVE_SPEED_DEFAULT
from odemis.acq.stream import UNDEFINED_ROI
from odemis.dataio import find_fittest_converter
from odemis.util import executeAsyncTask, dataio

ANCHOR_MAX_PIXELS = 512 ** 2  # max number of pixels for the anchor region


class MillingSettings(object):
    """
    Model class for milling settings

    """

    def __init__(self, name, current, horizontal_fov, roi,
                 pixel_size, duration, beam_angle=None,
                 dc_roi=UNDEFINED_ROI, dc_period=None, dc_dwell_time=None, dc_current=None):
        """
       Settings class for milling a rectangle.

       :param name: (str) name for the given milling setting
       :param current: (float) Probe current, as available in the presets, based on the preset naming. The ion beam current in A
       :param horizontal_fov:(float) width of the FoV when milling (to be set on the ion-beam settings), keep it fixed below 45um
       :param roi: (float, float, float, float) Region of interest, relative to the FoV (0→ 1)
       :param pixel_size:  float, float) Distance between the points, in m (Note: the orsay API expects a
        “probesize + overlap”, but we simplify it by just a pixelSize X/Y).
       :param beam_angle: (float) the angle between the stage and the beam column (for now this will be fixed, to ~10°), in rad
       :param duration: (float) total time to mill the region
       :param dc_roi: (float, float, float, float) Anchor region for the drift correction relative to the FoV (0→ 1).
        Use UNDEFINED_ROI if drift correction shouldn’t be applied
       :param dc_period: (float) time of milling before running the drift correction in s
       :param dc_dwell_time: (float) dwell time to be used during anchor region acquisition
       :param dc_current: (float) drift correction current
       :return:
       """
        self.horizontalFoV = model.FloatContinuous(horizontal_fov, unit="m", range=(5e-06, 45e-06))
        self.current = model.FloatContinuous(current, unit="A", range=(0.5e-12, 3000e-12))
        self.name = model.StringVA(name)
        self.roi = model.TupleContinuous(tuple(roi),
                                         range=((0, 0, 0, 0), (1, 1, 1, 1)),
                                         cls=(int, float))
        self.pixelSize = model.TupleContinuous(tuple(pixel_size), unit="m",
                                               range=((0, 0), (1e-3, 1e-3)),  # max 1 mm: arbitrary gigantic value
                                               cls=(int, float))
        # angles ot of this range would be dangerous for the hardware and probably not useful for the user
        if beam_angle is None:
            stage = model.getComponent(role="stage")
            stage_md = stage.getMetadata()
            try:
                self._ion_beam_angle = stage_md[model.MD_ION_BEAM_TO_SAMPLE_ANGLE]
            except KeyError:
                raise ValueError("Ion beam angle not defined in stage metadata")
        else:
            self.beamAngle = model.FloatContinuous(beam_angle, unit="rad", range=(math.radians(6), math.radians(10)))
        self.duration = model.FloatContinuous(duration, unit="s", range=(0.1, 100e3))
        # drift correction settings
        if dc_roi == UNDEFINED_ROI:
            if dc_period is None:
                dc_period = duration
            if dc_dwell_time is None:
                dc_dwell_time = 1e-09
            if dc_current is None:
                dc_current = 1e-12
        else:
            if dc_period is None:
                raise ValueError("dc_period has to be provided.")
            if dc_dwell_time is None:
                raise ValueError("dc_dwell_time has to be provided.")
            if dc_current is None:
                raise ValueError("dc_current has to be provided.")
        self.dcRoi = model.TupleContinuous(tuple(dc_roi),
                                           range=((0, 0, 0, 0), (1, 1, 1, 1)),
                                           cls=(int, float))
        self.dcPeriod = model.FloatContinuous(dc_period, unit="s", range=(0.1, 100e3))
        self.dcDwellTime = model.FloatContinuous(dc_dwell_time, unit="s", range=(1e-09, 1e-03))
        self.dcCurrent = model.FloatContinuous(dc_current, unit="A", range=(0.5e-12, 3000e-12))


def load_config(yaml_filename):
    """
    Load user input from yaml settings file.
    :param yaml_filename: (str) Filename path of user configuration file.
    :return: (list) List containing MillingSettings.
    """
    try:
        with open(yaml_filename, "r") as f:
            settings = yaml.safe_load(f)

    except yaml.YAMLError as exc:
        logging.error("Syntax error in milling settings yaml file: %s", exc)

    millings = []
    for ms in settings:
        milling_setting = MillingSettings(**ms)
        millings.append(milling_setting)

    return millings


# To handle the timeout error when the stage is not able to move to the desired position
# It logs the message and raises the MoveError exception
class MoveError(Exception):
    pass


class MillingRectangleTask(object):
    """
    This class represents a milling Task for milling rectangular regions on the sample.
    """

    def __init__(self, future: futures.Future, millings: list, sites: list, feature_post_status: str, acq_streams,
                 ebeam, sed, stage, aligner, log_path=None):
        """
        Constructor
        :param future: (ProgressiveFuture) the future that will be executing the task
        :param millings: (list of MillingSettings) Settings corresponding to each milling, to be milled in order
        :param sites: (list of Features) Each Feature to be milled
        :param feature_post_status: (str) value to set on the Feature at the end of a complete milling series
        :param acq_streams: type of acquisition streams to be used for milling
        :param ebeam: model component for the scanner
        :param sed: model component for the detector
        :param stage: model component for the stage
        :param aligner: model component for the aligner
        :param log_path: (str) path to the log anchor region acquisition used in drift correction
        """
        self._stage = stage
        self._scanner = ebeam
        self._detector = sed
        self._aligner = aligner

        self._future = future
        if future is not None:
            self._future.running_subf = model.InstantaneousFuture()
            self._future._task_lock = threading.Lock()

        self._settings = millings
        self._feature_status = feature_post_status  # site and feature means the same
        self.sites = sites
        self.streams = acq_streams

        self._log_path = log_path
        if log_path:
            filename = os.path.basename(self._log_path)
            if not filename:
                raise ValueError("Filename is not found on log path.")
            self._exporter = find_fittest_converter(filename)
            self._fn_bs, self._fn_ext = dataio.splitext(filename)
            self._log_dir = os.path.dirname(self._log_path)

        # internal class requirements between functions
        self._pass_duration = []  # one value per milling setting, duration(s) after which drift acquisition starts
        self._iterations = []  # one value per milling setting, numbers of scans of .roi
        stage_md = stage.getMetadata()
        try:
            self._ion_beam_angle = stage_md[model.MD_ION_BEAM_TO_SAMPLE_ANGLE]
        except KeyError:
            raise ValueError("Ion beam angle not defined in stage metadata")

        # Rough estimate of the stage movement speed, for estimating the extra
        # duration due to movements (copied from _tiledacq.py)
        self._move_speed = MOVE_SPEED_DEFAULT
        if model.hasVA(stage, "speed"):
            try:
                self._move_speed = (stage.speed.value["x"] + stage.speed.value["y"]) / 2
            except Exception as ex:
                logging.warning("Failed to read the stage speed: %s", ex)

        # estimate milling time from the list of milling settings-> millings
        time_estimate = 0
        drift_acq_time = 0

        for setting_nb, setting in enumerate(self._settings):
            time_estimate += setting.duration.value
            if setting.dcRoi.value != UNDEFINED_ROI:
                self._iterations.append(int(numpy.ceil(setting.duration.value / setting.dcPeriod.value)))
                self._pass_duration.append(setting.duration.value / self._iterations[setting_nb])

                drift_acq_time += self.estimate_drift_time(setting)
            else:
                self._iterations.append(1)
                self._pass_duration.append(setting.duration.value)

        self._milling_time_estimate = time_estimate + drift_acq_time  # time_estimate per feature

    def cancel(self, future: "Future") -> bool:
        """
        Canceler of acquisition task.
        :param future: the future that will be executing the task
        :return: True if it successfully cancelled (stopped) the future
        """
        logging.debug("Canceling milling procedure...")

        with future._task_lock:
            if future._task_state == FINISHED:
                return False
            future._task_state = CANCELLED
            future.running_subf.cancel()
            logging.debug("Milling procedure cancelled.")
        return True

    def estimate_milling_time(self, sites_done: int = 0, actual_time_per_site: float = None) -> float:
        """
        Estimates the milling time for the given feature.
        :param sites_done: number of sites already milled
        :param actual_time_per_site: actual milling time measured for a single site
        :return: (float > 0): the estimated time is in seconds
        """
        remaining_sites = (len(self.sites) - sites_done)
        if actual_time_per_site:
            milling_time = actual_time_per_site * remaining_sites
        else:
            milling_time = self._milling_time_estimate * remaining_sites

        # Time it takes to go from one site to the other
        prev_stage_pos_ref = None  # previous stage position used for reference in estimating time
        stage_time = 0
        for site in self.sites:
            stage_time += self.estimate_stage_movement_time(site, prev_stage_pos_ref)
            prev_stage_pos_ref = site.pos.value

        return milling_time + stage_time

    def estimate_stage_movement_time(self, site: CryoFeature, stage_pos_ref: tuple = None) -> float:
        """
        Estimation the time taken by the stage to move from current position to the location of given site
        :param site: (CryoFeature) consists X,Y,Z coordinates for stage location in m.
        :param stage_pos_ref: (tuple) reference position of the stage in m from where the stage moves.
        :return: (float) time to move the stage between two points in s.
        """
        # current position from x,y,z of stage position and eliminating rx,ry,rz
        if stage_pos_ref is None:
            stage_pos = self._stage.position.value
            current_pos = [stage_pos[an] for an in ("x", "y", "z")]
        else:
            current_pos = stage_pos_ref
        target_pos = site.pos.value  # list
        diff = [abs(target - current) for target, current in zip(target_pos, current_pos)]
        stage_time = math.sqrt(sum(d ** 2 for d in diff)) / self._move_speed

        return stage_time

    def estimate_drift_time(self, setting: MillingSettings) -> float:
        """
        Estimate the time taken to acquire drift correction images
        :param setting: (MillingSettings) milling settings
        return (float): estimated time to acquire 1 anchor area
        """
        if setting.dcRoi.value == UNDEFINED_ROI:
            return 0

        nb_anchor_scanning = math.ceil(setting.duration.value / setting.dcPeriod.value)
        drift_estimation = AnchoredEstimator(self._scanner, self._detector,
                                             setting.dcRoi.value, setting.dcDwellTime.value, ANCHOR_MAX_PIXELS,
                                             follow_drift=False)
        return drift_estimation.estimateAcquisitionTime() * nb_anchor_scanning

    def _move_to_site(self, site: CryoFeature):
        """
        Move the stage to the given site.
        :param site: (CryoFeature) The site to move to.
        :raises MoveError: if the stage failed to move to the given site.
        """
        target_pos = {"x": site.pos.value[0],
                      "y": site.pos.value[1],
                      "z": site.pos.value[2]}
        logging.debug("For feature %s moving the stage to %s m", site.name.value, target_pos)
        self._future.running_subf = self._stage.moveAbs(target_pos)
        stage_time = self.estimate_stage_movement_time(site)
        t = stage_time * 10 + 5  # adding extra margin
        try:
            self._future.running_subf.result(t)
        except TimeoutError:
            self._future.running_subf.cancel()
            raise MoveError(f"Failed to move the stage for feature {site.name.value} within {t} s")

        # The stage never *exactly* reaches the target position. => Store how
        # far it was away, according to the stage encoders. We could try to
        # compensate using beam shift. However, on the MIMAS stage, there is
        # currently too much imprecision that the encoders do not detect, so we
        # just log the information.
        actual_pos = self._stage.position.value
        diff_pos = [target_pos[a] - actual_pos[a] for a in ("x", "y", "z")]
        logging.debug("Stage reached %s, away from target by %s m", actual_pos, diff_pos)

        # Reset the beam shift to 0,0 in order to start from scratch with drift
        # compensation, to reduce the chances of reaching the limits of the shift.
        self._scanner.shift.value = 0, 0

    def _mill_all_settings(self, site: CryoFeature):
        """
        Iterates over all the milling settings for one site
        :param site: The site to mill.
        :raises MoveError: if the stage failed to move to the given site, or the beam shift failed to compensate drift.
        """
        self._move_to_site(site)

        for setting_nb, milling_settings in enumerate(self._settings):
            # Tilt stage angle
            stage_rx = self._ion_beam_angle - milling_settings.beamAngle.value
            logging.debug("Tilting stage to %f °", math.degrees(stage_rx))
            self._future.running_subf = self._stage.moveAbs({"rx": stage_rx})

            try:
                t = 15  # s
                self._future.running_subf.result(timeout=t)
            except TimeoutError:
                logging.debug("Failed to tilt the stage for feature %s within %s s", site.name.value, t)
                raise

            self._scanner.horizontalFoV.value = milling_settings.horizontalFoV.value

            # Make the blanker automatic (ie, disabled when acquiring)
            self._scanner.blanker.value = None

            # Initialize the drift corrector only if it is requested by the user
            if milling_settings.dcRoi.value != UNDEFINED_ROI:
                # Change the current to the drift correction current
                self._scanner.probeCurrent.value = milling_settings.dcCurrent.value
                drift_est = AnchoredEstimator(self._scanner, self._detector,
                                              milling_settings.dcRoi.value,
                                              milling_settings.dcDwellTime.value,
                                              max_pixels=ANCHOR_MAX_PIXELS, follow_drift=False)
                # acquire an image at the given location (RoI)
                drift_est.acquire()
                da = drift_est.raw[-1]
                pxs_anchor = da.metadata[model.MD_PIXEL_SIZE]
                if self._log_path:
                    fn = f"{self._fn_bs}-dc-{site.name.value}-{milling_settings.name.value}-0{self._fn_ext}"
                    self._exporter.export(os.path.join(self._log_dir, fn), [da])

            # Compute the probe size and overlap, based on the requested pixel size
            # For the probe size, use the pixel size in the largest dimension (typically X)
            # To adjust the pixel size in Y, compute the overlap so that probe_size * (1-overlap)) == pixelSize.
            # X is horizontal axis, Y is vertical axis
            pxs = milling_settings.pixelSize.value
            probe_size = max(pxs)
            overlap = [1 - (pxs[0] / probe_size), 1 - (pxs[1] / probe_size)]  # Always >= 0

            logging.debug(f"Milling setting: {milling_settings.name.value} for feature: {site.name.value}")

            self._scanner.probeCurrent.value = milling_settings.current.value

            if milling_settings.dcRoi.value == UNDEFINED_ROI:

                with self._future._task_lock:
                    if self._future._task_state == CANCELLED:
                        raise CancelledError()

                    self._future.running_subf = mill_rectangle(milling_settings.roi.value, self._scanner,
                                                               iteration=self._iterations[setting_nb],
                                                               duration=self._pass_duration[setting_nb],
                                                               probe_size=probe_size, overlap=overlap)
                try:
                    self._future.running_subf.result()
                    logging.debug(f"Milling {milling_settings.name.value} for feature {site.name.value} finished")
                except Exception as exp:
                    logging.error(
                        f"Milling setting {milling_settings.name.value} for feature {site.name.value} failed: {exp}")
                    raise

            else:
                for itr in range(self._iterations[setting_nb]):
                    with self._future._task_lock:
                        if self._future._task_state == CANCELLED:
                            raise CancelledError()

                        self._scanner.probeCurrent.value = milling_settings.current.value
                        self._future.running_subf = mill_rectangle(milling_settings.roi.value, self._scanner,
                                                                   iteration=1,
                                                                   duration=self._pass_duration[setting_nb],
                                                                   probe_size=probe_size, overlap=overlap)
                    try:
                        self._future.running_subf.result()
                        logging.debug(f"Milling {milling_settings.name.value} for feature {site.name.value} iteration {itr} finished")
                    except Exception as exp:
                        logging.error(
                            f"Milling {milling_settings.name.value} for feature {site.name.value} at iteration {itr} failed: {exp}")
                        raise

                    # Change the current to the drift correction current
                    self._scanner.probeCurrent.value = milling_settings.dcCurrent.value
                    # Acquire an image at the given location (RoI)
                    drift_est.acquire()

                    # Estimate the drift since the last correction
                    drift_est.estimate()
                    drift = (pxs_anchor[0] * drift_est.drift[0],
                             -(pxs_anchor[1] * drift_est.drift[1]))

                    # Move FIB to compensate drift
                    previous_shift = self._scanner.shift.value
                    shift = (drift[0] + previous_shift[0], drift[1] + previous_shift[1])  # m
                    self._scanner.shift.value = self._scanner.shift.clip(shift)
                    if self._scanner.shift.value != shift:  # check if it has been clipped
                        raise MoveError(f"Failed to set the beam shift to {shift} m, limited to {self._scanner.shift.value}")
                    logging.debug("Ion-beam shift in m: %s changed to: %s", previous_shift, self._scanner.shift.value)
                    if self._log_path:
                        fn = f"{self._fn_bs}-dc-{site.name.value}-{milling_settings.name.value}-{itr + 1}{self._fn_ext}"
                        self._exporter.export(os.path.join(self._log_dir, fn), drift_est.raw[-1])

            self._scanner.blanker.value = True
            with self._future._task_lock:
                if self._future._task_state == CANCELLED:
                    raise CancelledError()

    def _acquire_feature(self, site: CryoFeature):
        """
        Acquire the data for the given feature and add to the given site.
        :param site: (CryoFeature) The feature to acquire.
        """
        data = []
        try:
            # check if cancellation happened while the acquiring future is working
            with self._future._task_lock:
                if self._future._task_state == CANCELLED:
                    raise CancelledError()
                self._future.running_subf = acquire(self.streams)
            data, exp = self._future.running_subf.result()
            if exp:
                logging.error(
                    f"Acquisition for feature {site.name.value} partially failed: {exp}")
        except Exception as exp:
            logging.error(f"Acquisition for feature {site.name.value} failed: {exp}")

        # Check on the acquired data
        if not data:
            logging.warning("The acquired data array in stream %s for feature %s is empty", self.streams,
                            site.name)
        else:
            # Convert the data to StaticStreams, and add them to the feature
            new_streams = dataio.data_to_static_streams(data)
            site.streams.value.extend(new_streams)
            logging.info("The acquisition for stream %s for feature %s is done", self.streams, site.name)

    def run(self):
        """
        The main function of the task class, which will be called by the future asynchronously
        """
        self._future._task_state = RUNNING

        try:
            actual_time_per_site = None
            self._scanner.shift.value = (0, 0)  # reset drift correction to have more margin

            microscope = model.getMicroscope()
            posture_manager = MicroscopePostureManager(microscope)
            for site_idx, site in enumerate(self.sites):
                # Update progress on the milling sites left
                start_time = time.time()
                remaining_t = self.estimate_milling_time(site_idx, actual_time_per_site)
                self._future.set_end_time(time.time() + remaining_t)

                with self._future._task_lock:
                    if self._future._task_state == CANCELLED:
                        raise CancelledError()
                    logging.debug(f"Retracting the objective to set the imaging in FIB mode")
                    self._future.running_subf = posture_manager.cryoSwitchSamplePosition(MILLING)

                self._future.running_subf.result()

                # For one given site, move the stage at the site and
                # do milling with all milling settings for the given site
                # If timeout exception is raised during stage movement, continue milling the next site
                try:
                    self._mill_all_settings(site)
                except MoveError as ex:
                    logging.warning("Feature %s: %s. Skipping to next feature.", site.name.value, ex)
                    continue

                # Update the status of milling for the given site
                site.status.value = self._feature_status
                logging.debug(f"The milling of feature {site.name.value} completed")

                # Acquire the streams of the milled site
                self._acquire_feature(site)

                # Store the actual time during milling one site computed after moving the stage
                actual_time_per_site = time.time() - start_time

        except CancelledError:
            logging.debug("Stopping because milling was cancelled")
            raise
        except Exception:
            logging.exception("The milling failed")
            raise
        finally:
            # activate the blanker
            self._scanner.blanker.value = True
            self._future._task_state = FINISHED


def mill_features(millings: list, sites: list, feature_post_status, acq_streams, ebeam, sed, stage,
                  aligner, log_path=None) -> futures.Future:
    """
    Mill features on the sample.
    :param millings: (list of MillingSettings) Settings corresponding to each milling, to be milled in order
    :param sites: (list of Features) Each Feature to be milled
    :param feature_post_status: (str) value to set on the Feature at the end of a complete milling series
    :param acq_streams: type of acquisition streams to be used for milling
    :param ebeam: model component for the scanner
    :param sed: model component for the detector
    :param stage: model component for the stage
    :param aligner: model component for the objective
    :return: ProgressiveFuture
    """
    # Create a progressive future with running sub future
    future = model.ProgressiveFuture()
    # create acquisition task
    milling_task = MillingRectangleTask(future, millings, sites, feature_post_status, acq_streams, ebeam, sed, stage,
                                        aligner, log_path)
    # add the ability of cancelling the future during execution
    future.task_canceller = milling_task.cancel

    # set the progress of the future
    total_duration = milling_task.estimate_milling_time()
    future.set_end_time(time.time() + total_duration)

    # assign the acquisition task to the future
    executeAsyncTask(future, milling_task.run)

    return future


def estimate_milling_time(*args, **kwargs) -> float:
    """
    Estimate the duration of milling.
    :params: arguments are the same as mill_features() arguments.
    :return: (float > 0): estimated milling time in seconds
    """
    milling_task = MillingRectangleTask(None, *args, **kwargs)
    return milling_task.estimate_milling_time()


class MillingSettings2:
    """Represents milling settings for a single milling task"""

    def __init__(self, current: float, voltage: float, field_of_view: float, mode: str = "Serial", channel: str = "ion"):
        self.current = model.FloatContinuous(current, unit="A", range=(20e-12, 120e-9))  # TODO: migrate to float enum after testing
        self.voltage = model.FloatContinuous(voltage, unit="V", range=(0, 30e3))         # TODO: migrate to float enum after testing
        self.field_of_view = model.FloatContinuous(field_of_view, unit="m", range=(50e-06, 960e-06))
        self.mode = model.StringEnumerated(mode, choices=set(["Serial", "Parallel"]))
        self.channel = model.StringEnumerated(channel, choices=set(["ion"])) # TODO: add support for electron milling

    def to_json(self) -> dict:
        return {"current": self.current.value,
                "voltage": self.voltage.value,
                "field_of_view": self.field_of_view.value,
                "mode": self.mode.value,
                "channel": self.channel.value}

    @staticmethod
    def from_json(data: dict) -> "MillingSettings2":
        return MillingSettings2(current=data["current"],
                                voltage=data["voltage"],
                                field_of_view=data["field_of_view"],
                                mode=data.get("mode", "Serial"),
                                channel=data.get("channel", "ion"))

    def __repr__(self):
        return f"{self.to_json()}"

# TODO: migrate these to ABC, dataclass for Rectangle, Line, Circle, etc.
class MillingPatternParameters:
    """Represents milling pattern parameters"""

    def __init__(self, width: float, height: float, depth: float, rotation: float = 0.0, center = (0, 0), scan_direction: str = "TopToBottom", name: str = "Rectangle"):
        self.name = model.StringVA(name)
        self.width = model.FloatContinuous(width, unit="m", range=(1e-9, 900e-6))
        self.height = model.FloatContinuous(height, unit="m", range=(1e-9, 900e-6))
        self.depth = model.FloatContinuous(depth, unit="m", range=(1e-9, 100e-6))
        self.rotation = model.FloatContinuous(rotation, unit="rad", range=(0, 2 * math.pi))
        self.center = model.TupleContinuous(center, unit="m", range=((-1e3, -1e3), (1e3, 1e3)), cls=(int, float))
        self.scan_direction = model.StringEnumerated(scan_direction, choices=set(["TopToBottom", "BottomToTop", "LeftToRight", "RightToLeft"]))

    def to_json(self) -> dict:
        return {"name": self.name.value,
                "width": self.width.value,
                "height": self.height.value,
                "depth": self.depth.value,
                "rotation": self.rotation.value,
                "center_x": self.center.value[0],
                "center_y": self.center.value[1],
                "scan_direction": self.scan_direction.value}
    
    @staticmethod
    def from_json(data: dict):
        return MillingPatternParameters(width=data["width"], 
                                        height=data["height"], 
                                        depth=data["depth"], 
                                        rotation=data.get("rotation", 0), 
                                        center=(data.get("center_x", 0), data.get("center_y", 0)), 
                                        scan_direction=data.get("scan_direction", "TopToBottom"), 
                                        name=data.get("name", "Rectangle"))

    def __repr__(self):
        return f"{self.to_json()}"

class MillingTaskSettings:
    milling: MillingSettings2           # MillingSettings2
    patterns: list                      # list[MillingPatternParameters]
    drift_correction: dict              # NotYetImplemented

    def __init__(self, milling: dict, patterns: list):
        self.milling = milling
        self.patterns = patterns

    def to_json(self) -> dict:
        return {"milling": self.milling.to_json(), "patterns": [pattern.to_json() for pattern in self.patterns]}

    @staticmethod
    def from_json(data: dict):
        return MillingTaskSettings(milling=MillingSettings2(**data["milling"]), 
                                   patterns=[MillingPatternParameters.from_json(pattern) for pattern in data["patterns"]])

    def __repr__(self):
        return f"{self.to_json()}"

class MillingTask:
    """This class represents a standard Milling Task."""

    def __init__(self, future: futures.Future, settings: MillingTaskSettings):
        """
        :param future: (ProgressiveFuture) the future that will be executing the task
        :param settings: (MillingTaskSettings) the settings representing the milling task
        """

        self.microscope = model.getComponent(role="fibsem")
        self.settings = settings

        self._future = future
        if future is not None:
            self._future.running_subf = model.InstantaneousFuture()
            self._future._task_lock = threading.Lock()

    def cancel(self, future: "Future") -> bool:
        """
        Canceler of acquisition task.
        :param future: the future that will be executing the task
        :return: True if it successfully cancelled (stopped) the future
        """
        logging.debug("Canceling milling procedure...")

        with future._task_lock:
            if future._task_state == FINISHED:
                return False
            future._task_state = CANCELLED
            future.running_subf.cancel()
            self.microscope.stop_milling()
            logging.debug("Milling procedure cancelled.")
        return True

    def estimate_milling_time(self) -> float:
        """
        Estimates the milling time for the given patterns.
        :return: (float > 0): the estimated time is in seconds
        """
        return self.microscope.estimate_milling_time()

    def run_milling(self, settings: MillingTaskSettings):
        """Run the milling task with the given settings. ThermoFisher implementation"""
        microscope = self.microscope

        # get the milling settings
        milling_current = settings.milling.current.value
        milling_voltage = settings.milling.voltage.value
        milling_fov = settings.milling.field_of_view.value
        milling_channel = settings.milling.channel.value
        milling_mode = settings.milling.mode.value

        # get initial imaging settings
        imaging_current = microscope.get_beam_current(milling_channel)
        imaging_voltage = microscope.get_high_voltage(milling_channel)
        imaging_fov = microscope.get_field_of_view(milling_channel)

        # error management
        ce = False
        try:

            # set the milling state
            microscope.clear_patterns()
            microscope.set_default_patterning_beam_type(milling_channel)
            microscope.set_high_voltage(milling_voltage, milling_channel)
            microscope.set_beam_current(milling_current, milling_channel)
            microscope.set_field_of_view(milling_fov, milling_channel)
            microscope.set_patterning_mode(milling_mode)

            # draw milling patterns to screen
            for pattern in settings.patterns:
                microscope.create_rectangle(pattern.to_json())

            # estimate the milling time
            estimated_time = microscope.estimate_milling_time()
            self._future.set_end_time(time.time() + estimated_time)

            # start patterning (async)
            microscope.start_milling()

            # wait for milling to finish
            elapsed_time = 0
            wait_time = 5
            while microscope.get_patterning_state() == "Running":

                with self._future._task_lock:
                    if self._future.cancelled() == CANCELLED:
                        raise CancelledError()

                logging.info(f"Milling in progress... {elapsed_time} / {estimated_time}")
                time.sleep(wait_time)
                elapsed_time += wait_time

        except CancelledError as ce:
            logging.info(f"Cancelled milling: {ce}")
        except Exception as e:
            logging.error(f"Error while milling: {e}")
        finally:
            # restore imaging state
            microscope.set_beam_current(imaging_current, milling_channel)
            microscope.set_high_voltage(imaging_voltage, milling_channel)
            microscope.set_field_of_view(imaging_fov, milling_channel)
            microscope.set_active_view(2) # ion beam TODO: use milling_channel
            microscope.clear_patterns()

            if ce:
                raise ce # future excepts error to be raised
        return

    def run(self):
        """
        The main function of the task class, which will be called by the future asynchronously
        """
        self._future._task_state = RUNNING

        try:

            with self._future._task_lock:
                if self._future._task_state == CANCELLED:
                    raise CancelledError()

            self.run_milling(self.settings)
            logging.debug("The milling completed")

        except CancelledError:
            logging.debug("Stopping because milling was cancelled")
            raise
        except Exception:
            logging.exception("The milling failed")
            raise
        finally:
            self._future._task_state = FINISHED


def mill_patterns(settings: MillingTaskSettings) -> futures.Future:
    """
    Run Mill patterns.
    :param settings: (MillingTaskSettings) Settings for the milling task
    :return: ProgressiveFuture
    """
    # Create a progressive future with running sub future
    future = model.ProgressiveFuture()
    # create acquisition task
    milling_task = MillingTask(future, settings)
    # add the ability of cancelling the future during execution
    future.task_canceller = milling_task.cancel

    # set the progress of the future (TODO: fix dummy time estimate)
    future.set_end_time(time.time() + 10 * len(settings.patterns))

    # assign the acquisition task to the future
    executeAsyncTask(future, milling_task.run)

    return future
