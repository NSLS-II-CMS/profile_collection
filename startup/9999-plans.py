print(f"Loading {__file__!r} ...")

import time
import bluesky.preprocessors as bpp
from bluesky.utils import short_uid
from functools import partial

try:
    DETS = get_beamline().detector + [core_laser, laserx, lasery, smy, smx, sth, schi]
except Exception:
    DETS = get_beamline().detector + [smy, smx, sth, schi]

sample_pta = Sample('test')

@bpp.finalize_decorator(final_plan=shutter_off)
def expose(detectors, exposure_time=None, extra=None, verbosity=3, md=None):
    md = dict(md or {})
    md.setdefault("measure_type", "expose")

    # Set exposure time
    if exposure_time is not None:
        exposure_time = abs(exposure_time)

        for detector in detectors:
            if (
                hasattr(detector, "cam")
                and hasattr(detector, "setExposureTime")
                and exposure_time != detector.cam.acquire_time.get()
            ):
                yield from detector.setExposureTime(exposure_time, verbosity=verbosity)

    yield from shutter_on()

    md["plan_header_override"] = md["measure_type"]
    start_time = time.time()

    # md_current = yield from self.get_md()
    md["beam_int_bim3"] = yield from beam.bim3.flux(verbosity=0)
    md["beam_int_bim4"] = yield from beam.bim4.flux(verbosity=0)
    md["beam_int_bim5"] = yield from beam.bim5.flux(verbosity=0)
    # md['trigger_time'] = self.clock()
    # md.update(md_current)

    uids = yield from count(detectors, md=md)


def measure_single(
    sample,
    detectors=None,
    exposure_time=None,
    measure_type="measure",
    verbosity=3,
    handlefile=True,
    md=None,
):
    """Measure data by triggering the area detectors.

    Parameters
    ----------
    exposure_time : float
        How long to collect data
    extra : string, optional
        Extra information about this particular measurement (which is typically
        included in the savename/filename).
    """
    detectors = detectors if detectors is not None else DETS
    md = dict(md or {})
    if exposure_time is not None:
        sample.set_attribute("exposure_time", exposure_time)
    # else:
    # exposure_time = sample.get_attribute('exposure_time')
    extra = md.get('extra', None)
    savename = sample.get_savename(savename_extra=extra)

    if verbosity >= 2 and (get_beamline().current_mode != "measurement"):
        print("WARNING: Beamline is not in measurement mode (mode is '{}')".format(get_beamline().current_mode))

    if verbosity >= 1 and len(get_beamline().detector) < 1:
        raise ValueError("ERROR: No detectors defined in detectors")

    md_current = yield from sample.get_md()
    new_md = yield from sample.get_measurement_md()
    md_current.update(new_md)
    md_current["sample_savename"] = savename
    md_current["measure_type"] = measure_type
    # md_current['filename'] = '{:s}_{:04d}.tiff'.format(savename, md_current['detector_sequence_ID'])
    # md_current['filename'] = '{:s}_{:04d}.tiff'.format(savename, RE.md['scan_id'])
    md_current["filename"] = "{:s}_{:06d}".format(savename, RE.md.get("scan_id", -1))
    md_current.update(md)

    yield from expose(detectors, exposure_time, extra=extra, verbosity=verbosity, md=md_current)
    # sample.expose(exposure_time, extra=extra, verbosity=verbosity, **md)

    # This is the symlinking code.
    # We plan to remove this once kafka based linking is finished.
    if handlefile:
        for detector in detectors:
            sample.handle_file(detector, verbosity=verbosity, **md_current)

    sample.md["measurement_ID"] += 1


def tiling(detectors, inner_plan, tiling_type=None, md=None):
    """
    A helper function that applies tiling to a plan.
    Tiling is used to fill gaps between the detector chips.
    
    There are 3 different tiling modes: None, xgaps, and xygaps.
    
    xgaps mode sets the detector positions to the upper position, calls the plan,
    and then sets the detector position to the lower position and calls the plan again.

    xygaps mode executes the plan for four different detector postions.

    The motors that are used to tiling will be added on to the list of detectors
    when passed to the inner plan.

    Parameters
    ----------
    detectors: list
        The list of detectors used in the plan.
        Detectors in this list that we know how to tiled, will get tiled.
    inner_plan: callable
        expected signature:

           def plan(detectors, md):
               '''
               Parameters
               ----------
               detectors: list
                   The list of detectors to use in the plan.
               md: dict
                   Plan metadata.
               '''
               ...
    tiling_type: str, None
        has one of the following values: None, 'xgaps', 'xygaps'
    md: dict
        Plan metadata.

    Returns
    -------
    list
        accumulated results of inner_plan

    """ 
    md = dict(md or {})
    md.setdefault('tile_id', short_uid('tile_id')) 
    GAP_SIZE = 5.16

    offsets = {
        "lower": {"saxs_x": 0, "saxs_y": 0, "waxs_x": 0, "waxs_y": 0},
        "upper": {"saxs_x": 0, "saxs_y": GAP_SIZE, "waxs_x": 0, "waxs_y": GAP_SIZE},
        "lower_left": {"saxs_x": 0, "saxs_y": 0, "waxs_x": 0, "waxs_y": 0},
        "upper_left": {"saxs_x": 0, "saxs_y": GAP_SIZE, "waxs_x": 0, "waxs_y": GAP_SIZE},
        "lower_right": {"saxs_x": GAP_SIZE, "saxs_y": 0, "waxs_x": -GAP_SIZE, "waxs_y": 0},
        "upper_right": {"saxs_x": GAP_SIZE, "saxs_y": GAP_SIZE, "waxs_x": -GAP_SIZE, "waxs_y": GAP_SIZE},
        "default": {"saxs_x": 0, "saxs_y": 0, "waxs_x": 0, "waxs_y": 0},
    }
    
    tile_types = {
        "xygaps": ["lower_left", "upper_left", "lower_right", "upper_right"],
        "ygaps": ["upper", "lower"],
        None: ["default"],
    }
    extra = md.get('extra', None)
    extras = {'lower': "pos1" if extra is None else f"{extra}_pos1",
              'upper': "pos2" if extra is None else f"{extra}_pos2",
              'lower_left': "pos1" if extra is None else f"{extra}_pos1",
              'upper_left': "pos2" if extra is None else f"{extra}_pos2",
              'lower_right': "pos3" if extra is None else f"{extra}_pos3",
              'upper_right': "pos4" if extra is None else f"{extra}_pos4",
              'default': extra
              }

    motors = []
    if pilatus2M in detectors:
        motors.extend([SAXSx, SAXSy])
    if pilatus800 in detectors:
        motors.extend([WAXSx, WAXSy])
    
    @bpp.reset_positions_decorator(motors)
    def tiling_wrapper():
        ret = []
        if pilatus2M in detectors:
            SAXSy_original = yield from bps.rd(SAXSy)
            SAXSx_original = yield from bps.rd(SAXSx)
        if pilatus800 in detectors:
            WAXSy_original = yield from bps.rd(WAXSy)
            WAXSx_original = yield from bps.rd(WAXSx)

        for tile in tile_types[tiling_type]:
            if pilatus2M in detectors:
                yield from bps.mv(SAXSx, SAXSx_original + offsets[tile]['saxs_x'])
                yield from bps.mv(SAXSy, SAXSy_original + offsets[tile]['saxs_y'])
            if pilatus800 in detectors:
                yield from bps.mv(WAXSx, WAXSx_original + offsets[tile]['waxs_x'])
                yield from bps.mv(WAXSy, WAXSy_original + offsets[tile]['waxs_y'])
            val = yield from inner_plan(detectors + motors, 
                                        md={**md, 'extra':extra, 'detector_position': tile})
            ret.append(val)
        return ret
    return (yield from tiling_wrapper())


def measure(
    sample,
    detectors=None,
    exposure_time=None,
    extra=None,
    measure_type="measure",
    verbosity=3,
    tiling_type=None,
    md=None,
):
    """Measure data by triggering the area detectors.

    Parameters
    ----------
    sample: Sample
        The sample object to measure.
    detectors: interable
        The list of detectors.
    exposure_time : float
        How long to collect data
    extra : string, optional
        Extra information about this particular measurement (which is typically
        included in the savename/filename).
    tiling : string
        Controls the detector tiling mode.
        None : regular measurement (single detector position)
        'ygaps' : try to cover the vertical gaps in the Pilatus detector
    """
    detectors = detectors if detectors is not None else DETS
    md = dict(md or {})

    yield from tiling(
            detectors,
            partial(measure_single,
                sample,
                exposure_time=exposure_time,
                measure_type=measure_type,
                verbosity=verbosity),
            tiling_type=tiling_type,
            md=md
            )
        

def sam_measure(
    sample=None,
    detectors=None,
    exposure_time=None,
    extra=None,
    measure_type="measure",
    verbosity=3,
    tiling=None,
    md=None,
):  
    sample = sample or sample_pta
    detectors = detectors if detectors is not None else DETS
    md = dict(md or {})
    if tiling is not None:
        raise ValueError("Parameter 'tiling' must be None. Other values are not supported yet.")

    yield from measure(
        sample,
        detectors=detectors,
        exposure_time=exposure_time,
        extra=extra,
        measure_type=measure_type,
        verbosity=verbosity,
        tiling=tiling,
        md=md,
    )
