import time
import bluesky.preprocessors as bpp

DETS = get_beamline().detector + [core_laser, laser, laserx, lasery, smy, smx, sth, schi]


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
    detectors=DETS,
    exposure_time=None,
    extra=None,
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
    md = dict(md or {})
    if exposure_time is not None:
        sample.set_attribute("exposure_time", exposure_time)
    # else:
    # exposure_time = sample.get_attribute('exposure_time')

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
            sample.handle_file(detector, extra=extra, verbosity=verbosity, md=md)

    sample.md["measurement_ID"] += 1


def measure(
    sample,
    detectors=DETS,
    exposure_time=None,
    extra=None,
    measure_type="measure",
    verbosity=3,
    tiling=None,
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
    GAP_SIZE = 5.16

    md = dict(md or {})

    offsets = {
        "lower": {"saxs_x": 0, "saxs_y": 0, "waxs_x": 0, "waxs_y": 0},
        "upper": {"saxs_x": 0, "saxs_y": GAP_SIZE, "waxs_x": 0, "waxs_y": GAP_SIZE},
        "lower_left": {"saxs_x": 0, "saxs_y": 0, "waxs_x": 0, "waxs_y": 0},
        "upper_left": {"saxs_x": 0, "saxs_y": GAP_SIZE, "waxs_x": 0, "waxs_y": GAP_SIZE},
        "lower_right": {"saxs_x": GAP_SIZE, "saxs_y": 0, "waxs_x": -GAP_SIZE, "waxs_y": 0},
        "upper_right": {"saxs_x": GAP_SIZE, "saxs_y": GAP_SIZE, "waxs_x": -GAP_SIZE, "waxs_y": GAP_SIZE},
        "default": {"saxs_x": 0, "saxs_y": 0, "waxs_x": 0, "waxs_y": 0},
    }

    positions = {
        "xygaps": ["lower_left", "upper_left", "lower_right", "upper_right"],
        "ygaps": ["upper", "lower"],
        None: ["default"],
    }

    extras = {'lower': "pos1" if extra is None else f"{extra}_pos1",
              'upper': "pos2" if extra is None else f"{extra}_pos2",
              'lower_left': "pos1" if extra is None else f"{extra}_pos1",
              'upper_left': "pos2" if extra is None else f"{extra}_pos2",
              'lower_right': "pos3" if extra is None else f"{extra}_pos3",
              'upper_right': "pos4" if extra is None else f"{extra}_pos4",
              }

    # TODO: Maybe this should raise if 2M and 800 are not in detectors, and tiling is not None.

    for position in positions[tiling]:

        SAXSy_original = yield from bps.rd(SAXSy)
        SAXSx_original = yield from bps.rd(SAXSx)
        WAXSy_original = yield from bps.rd(WAXSy)
        WAXSx_original = yield from bps.rd(WAXSx)

        if pilatus2M in detectors:
            yield from bps.mv(SAXSx, SAXSx_original + offsets[position]['saxs_x'])
            yield from bps.mv(SAXSy, SAXSy_original + offsets[position]['saxs_y'])
        if pilatus800 in detectors:
            yield from bps.mv(WAXSx, WAXSx_original + offsets[position]['waxs_x'])
            yield from bps.mv(WAXSy, WAXSy_original + offsets[position]['waxs_y'])

        md["detector_position"] = position
        extra_current = extras[position]

        yield from measure_single(
            sample,
            detectors=detectors,
            exposure_time=exposure_time,
            extra=extra_current,
            measure_type=measure_type,
            verbosity=verbosity,
            md=md,
        )
    
    if pilatus2M in detectors:
        yield from bps.mv(SAXSx, SAXSx_original)
        yield from bps.mv(SAXSy, SAXSy_original)
    if pilatus800 in detectors:
        yield from bps.mv(WAXSx, SAXSx_original)
        yield from bps.mv(WAXSy, SAXSy_original)


def sam_measure(
    sample=sample_pta,
    detectors=DETS,
    exposure_time=None,
    extra=None,
    measure_type="measure",
    verbosity=3,
    tiling=None,
    md=None,
):
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
