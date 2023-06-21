def sam_measure(
    exposure_time=None, extra=None, measure_type="measure", verbosity=3, tiling=None, stitchback=False, **md
):

    if tiling is not None:
        raise ValueError("Parameter 'tiling' must be None. Other values are not supported yet.")

    yield from sam.measure(
        exposure_time=exposure_time,
        extra=extra,
        measure_type=measure_type,
        verbosity=verbosity,
        tiling=tiling,
        stitchback=stitchback,
        **md,
    )


def measure_single(sample, exposure_time=None, extra=None, measure_type="measure", verbosity=3, **md):
    """Measure data by triggering the area detectors.

    Parameters
    ----------
    exposure_time : float
        How long to collect data
    extra : string, optional
        Extra information about this particular measurement (which is typically
        included in the savename/filename).
    """
    if exposure_time is not None:
        sample.set_attribute("exposure_time", exposure_time)
    # else:
    # exposure_time = sample.get_attribute('exposure_time')

    savename = sample.get_savename(savename_extra=extra)

    if verbosity >= 2 and (get_beamline().current_mode != "measurement"):
        print(
            "WARNING: Beamline is not in measurement mode (mode is '{}')".format(get_beamline().current_mode)
        )

    if verbosity >= 1 and len(get_beamline().detector) < 1:
        raise ValueError("ERROR: No detectors defined in cms.detector")

    md_current = yield from sample.get_md()
    new_md = yield from sample.get_measurement_md()
    md_current.update(new_md)
    md_current["sample_savename"] = savename
    md_current["measure_type"] = measure_type
    # md_current['filename'] = '{:s}_{:04d}.tiff'.format(savename, md_current['detector_sequence_ID'])
    # md_current['filename'] = '{:s}_{:04d}.tiff'.format(savename, RE.md['scan_id'])
    md_current["filename"] = "{:s}_{:06d}".format(savename, RE.md.get("scan_id", -1))
    md_current.update(md)

    yield from sample.expose(exposure_time, extra=extra, verbosity=verbosity, **md_current)
    # sample.expose(exposure_time, extra=extra, verbosity=verbosity, **md)

    sample.md["measurement_ID"] += 1



def measure(
    sample, exposure_time=None, extra=None, measure_type="measure", verbosity=3, tiling=None, stitchback=False, **md
):
        """Measure data by triggering the area detectors.

        Parameters
        ----------
        sample: Sample
            The sample object to measure.
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

        if tiling == "xygaps":
            SAXSy_o = SAXSy.user_readback.value
            SAXSx_o = SAXSx.user_readback.value
            WAXSy_o = WAXSy.user_readback.value
            WAXSx_o = WAXSx.user_readback.value
            # MAXSy_o = MAXSy.user_readback.value

            extra_current = "pos1" if extra is None else "{}_pos1".format(extra)
            md["detector_position"] = "lower_left"
            yield from measure_single(
                sample,
                exposure_time=exposure_time,
                extra=extra_current,
                measure_type=measure_type,
                verbosity=verbosity,
                stitchback=True,
                **md,
            )

            # pos2
            if [pilatus2M] in cms.detector:
                SAXSy.move(SAXSy_o + 5.16)
            if [pilatus800] in cms.detector:
                WAXSy.move(WAXSy_o + 5.16)
            # if [pilatus300] in cms.detector:
            # MAXSy.move(MAXSy_o + 5.16)

            extra_current = "pos2" if extra is None else "{}_pos2".format(extra)
            md["detector_position"] = "upper_left"
            yield from measure_single(
                sample, exposure_time=exposure_time, extra=extra_current, verbosity=verbosity, stitchback=True, **md
            )

            # pos4
            if pilatus2M in cms.detector:
                SAXSx.move(SAXSx_o + 5.16)
                SAXSy.move(SAXSy.o + 5.16)
            if pilatus800 in cms.detector:
                WAXSx.move(WAXSx_o - 5.16)
                WAXSy.move(WAXSy_o + 5.16)
            extra_current = "pos4" if extra is None else "{}_pos4".format(extra)
            md["detector_position"] = "upper_right"
            yield from measure_single(
                sample, exposure_time=exposure_time, extra=extra_current, verbosity=verbosity, stitchback=True, **md
            )

            # pos3
            if pilatus2M in cms.detector:
                SAXSx.move(SAXSx_o + 5.16)
                SAXSy.move(SAXSy_o)
            if pilatus800 in cms.detector:
                WAXSx.move(WAXSx_o - 5.16)
                WAXSy.move(WAXSy_o)

            extra_current = "pos3" if extra is None else "{}_pos3".format(extra)
            md["detector_position"] = "lower_right"
            yield from measure_single(
                sample, exposure_time=exposure_time, extra=extra_current, verbosity=verbosity, stitchback=True, **md
            )

            if WAXSx.user_readback.value != WAXSx_o:
                WAXSx.move(WAXSx_o)
            if WAXSy.user_readback.value != WAXSy_o:
                WAXSy.move(WAXSy_o)

            if SAXSx.user_readback.value != SAXSx_o:
                SAXSx.move(SAXSx_o)
            if SAXSy.user_readback.value != SAXSy_o:
                SAXSy.move(SAXSy_o)

        elif tiling == "ygaps":
            SAXSy_o = SAXSy.user_readback.value
            SAXSx_o = SAXSx.user_readback.value
            WAXSy_o = WAXSy.user_readback.value
            WAXSx_o = WAXSx.user_readback.value
            ##MAXSy_o = MAXSy.user_readback.value

            extra_current = "pos1" if extra is None else "{}_pos1".format(extra)
            md["detector_position"] = "lower"
            yield from measure_single(
                sample,
                exposure_time=exposure_time,
                extra=extra_current,
                measure_type=measure_type,
                verbosity=verbosity,
                stitchback=True,
                **md,
            )

            if pilatus2M in cms.detector:
                SAXSy.move(SAXSy_o + 5.16)
            if pilatus800 in cms.detector:
                WAXSy.move(WAXSy_o + 5.16)
            # if pilatus300 in cms.detector:
            # MAXSy.move(MAXSy_o + 5.16)

            # extra x movement is needed for pilatus2M.
            extra_current = "pos2" if extra is None else "{}_pos2".format(extra)
            md["detector_position"] = "upper"
            yield from measure_single(
                sample,
                exposure_time=exposure_time,
                extra=extra_current,
                measure_type=measure_type,
                verbosity=verbosity,
                stitchback=True,
                **md,
            )

            if SAXSy.user_readback.value != SAXSy_o:
                SAXSy.move(SAXSy_o)
            if WAXSy.user_readback.value != WAXSy_o:
                WAXSy.move(WAXSy_o)
            # if MAXSy.user_readback.value != MAXSy_o:
            # MAXSy.move(MAXSy_o)

        else:
            # Just do a normal measurement
            yield from measure_single(
                sample, exposure_time=exposure_time, extra=extra, measure_type=measure_type, verbosity=verbosity, **md
            )    
