#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

#######
# v3 -->  v4:  coord based on KY_coord_form.py

# applied on sample : 16(7th) , 7(8th), 6(9th), 11(10th), 12(11th), 1(12th, LonC), 5(13th), 3(14th, LonC)

#######


################################################################################
#  Short-term settings (specific to a particular user/experiment) can
# be placed in this file. You may instead wish to make a copy of this file in
# the user's data directory, and use that as a working copy.
################################################################################


from ophyd import EpicsSignal, Device, Component as Cpt
from bluesky.suspenders import SuspendFloor, SuspendCeil
from bluesky.preprocessors import stage_decorator

import sympy as sym

if True:
    # Define suspenders to hold data collection if x-ray
    # beam is not available.

    ring_current = EpicsSignal("SR:OPS-BI{DCCT:1}I:Real-I")
    sus = SuspendFloor(ring_current, 100, resume_thresh=400, sleep=600)
    RE.install_suspender(sus)

#    absorber_pos = EpicsSignal( 'XF:11BMB-ES{SM:1-Ax:ArmR}Mtr.RBV')
#    sus_abs_low = SuspendFloor(absorber_pos, -56, resume_thresh=-55)
#    sus_abs_hi = SuspendCeil(absorber_pos, -54, resume_thresh=-55)
#    RE.install_suspender(sus_abs_low)
#    RE.install_suspender(sus_abs_hi)

# RE.clear_suspenders()


if False:
    # The following shortcuts can be used for unit conversions. For instance,
    # for a motor operating in 'mm' units, one could instead do:
    #     sam.xr( 10*um )
    # To move it by 10 micrometers. HOWEVER, one must be careful if using
    # these conversion parameters, since they make implicit assumptions.
    # For instance, they assume linear axes are all using 'mm' units. Conversely,
    # you will not receive an error if you try to use 'um' for a rotation axis!
    m = 1e3
    cm = 10.0
    mm = 1.0
    um = 1e-3
    nm = 1e-6

    inch = 25.4
    pixel = 0.172  # Pilatus

    deg = 1.0
    rad = np.degrees(1.0)
    mrad = np.degrees(1e-3)
    urad = np.degrees(1e-6)


INTENSITY_EXPECTED_050 = 18800.0
INTENSITY_EXPECTED_025 = INTENSITY_EXPECTED_050 * 0.5


def get_default_stage():
    return stg

class Laser(Device):
    pulse_width = Cpt(EpicsSignal, 'Width-SP')
    cmd = Cpt(EpicsSignal, 'Cmd')
    turn_on = Cpt(EpicsSignal, 'Out-Sel')

laser = Laser('biome:{Trigger}', name="laser")

class SampleTSAXS(SampleTSAXS_Generic):
    def __init__(self, name, base=None, **md):
        super().__init__(name=name, base=base, **md)
        self.naming_scheme = ["name", "extra", "temperature", "exposure_time"]

        self.md["exposure_time"] = 30.0


class SampleGISAXS(SampleGISAXS_Generic):
    def __init__(self, name, base=None, **md):
        super().__init__(name=name, base=base, **md)
        self.naming_scheme = ["name", "extra", "th", "exposure_time"]


# class Sample(SampleTSAXS):
class Sample(SampleGISAXS):
    def __init__(self, name, base=None, **md):
        super().__init__(name=name, base=base, **md)

        # self.naming_scheme = ['name', 'extra', 'clock', 'temperature', 'th', 'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'temperature', 'th', 'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'th', 'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'y', 'th', 'clock', 'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'opos', 'lpos', 'x', 'th',  'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'clock', 'localT', 'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'x', 'yy',  'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'SamX', 'SamY',  'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'x', 'yy', 'SamX', 'SamY', 'exposure_time']
        # self.naming_scheme = ['name', 'extra', 'x', 'yy', 'h1', 'h2', 'exposure_time']
        self.naming_scheme = ["name", "extra", "Tc", "clock", "x", "th", "exposure_time"]
        self.name_o = self.name
        self.md["exposure_time"] = 5.0
        # self.incident_angles_default = [0.08, 0.10, 0.12, 0.15, 0.20]

        #       self.anneal_time = int(self.name.split('_')[-2].split('anneal')[-1])
        # self.anneal_time = int(self.name.split('anneal')[-1].split('_')[0])
        # self.preanneal_time = int(self.name.split('pre')[-1].split('_')[0])

        self._positional_axis = ["x", "y"]

        self.smxPos = [121.5, 142.5, 162.5]

        # self._axes["x"].origin = 130
        # self._axes["y"].origin = 16.8  # smy stage should be set with the limit [-5.5, -5]
        # self._axes["th"].origin = 0

        self._axes["x"].origin = 0
        self._axes["y"].origin = 39.26  # smy stage should be set with the limit [-5.5, -5]
        self._axes["th"].origin = 0.93

        #default position for laser position on the edge of the sample (5mm offset)
        #smx = -1.5, laserx = 0
        #smx and laserx should move simultaneously. 

    # def get_attribute(self, attribute):
    #     '''Return the value of the requested md.'''
    #     #return the value of temperature by IR laser where the
    #     if attribute == 'localT':
    #         return 'localT{:.1f}'.format(self.x2localT())
    #     elif attribute == 'clockT':
    #         return 'clock{:.2f}'.format(self.clock())
    #     else:
    #         return super().get_attribute(attribute)

    # #convert temperature to x position to tempreature in the gradient created by laser
    # def localT2x(self, temperature, t_range=[25, 500], xo=0, length=10):
    #     #xo is located at the HOT temperature (laser spot)
    #     #(x_pos-xo)/(-length) = (temperature-t_range[1])/(t_range[0]-t_range[1])
    #     #x_pos = (temperature-t_range[1])/(t_range[0]-t_range[1])*(-length) + xo

    #     x_pos = (492 - temperature)/24.7 # P_frac = 0.45
    #     return -x_pos

    # #convert x position to tempreature in the gradient created by laser
    # def x2localT(self, t_range=[25, 500], xo=0, length=10):
    #     #(x_pos-xo)/(-length) = (temperature-t_range[1])/(t_range[0]-t_range[1])
    #     x_pos = self.xpos(verbosity=0)
    #     #temperature = (x_pos-xo)/(-length)*(t_range[0]-t_range[1]) + t_range[1]

    #     T_est = 492 - abs(x_pos)*24.7 # P_frac = 0.45

    #     return T_est

    def setFlow(self, channel, voltage=0):
        # device = 'A1'
        ioL.set(AO[channel], 0)
        time.sleep(1)
        ioL.set(AO[channel], voltage)

    def setDryFlow(self, voltage=None):
        if voltage == None or voltage > 5 or voltage < 0:
            print("Input voltage betwee 0 and 5V")
        self.setFlow(1, voltage=voltage)

    def _set_axes_definitions(self):
        """Internal function which defines the axes for this stage. This is kept
        as a separate function so that it can be over-ridden easily."""

        # The _axes_definitions array holds a list of dicts, each defining an axis
        super()._set_axes_definitions()

        # self._axes_definitions.append ( {'name': 'y',
        #                     'motor': smy2,
        #                     'enabled': True,
        #                     'scaling': +1.0,
        #                     'units': 'mm',
        #                     'hint': 'positive moves stage up',
        #                     },
        #                     {'name': 'x',
        #                     'motor': ,
        #                     'enabled': True,
        #                     'scaling': +1.0,
        #                     'units': 'mm',
        #                     'hint': 'positive moves stage up',
        #                     } )

    def get_attribute(self, attribute):
        """Return the value of the requested md."""
        if attribute == "Tc":
            # return 'SamX{:.2f}'.format(-1*self.yypos(verbosity=0))     # Nov 2020
            return "Tc{:.2f}".format(pta.getTemperature(verbosity=0))  # Feb 2021
        else:
            return super().get_attribute(attribute)

    def get_md(self, prefix="sample_", include_marks=True, **md):
        """Returns a dictionary of the current metadata.
        The 'prefix' argument is prepended to all the md keys, which allows the
        metadata to be grouped with other metadata in a clear way. (Especially,
        to make it explicit that this metadata came from the sample.)"""

        # Update internal md
        # self.md['key'] = value

        yield from bps.null()  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        md_return = self.md.copy()
        md_return["name"] = self.name

        if include_marks:
            for label, positions in self._marks.items():
                md_return["mark_" + label] = positions

        # Add md that varies over time
        md_return["clock"] = self.clock()

        for axis_name, axis in self._axes.items():
            md_return[axis_name] = axis.get_position(verbosity=0)
            md_return["motor_" + axis_name] = axis.get_motor_position(verbosity=0)

        md_return["savename"] = self.get_savename()  # This should be over-ridden by 'measure'

        # Include the user-specified metadata
        md_return.update(md)

        # Add an optional prefix
        if prefix is not None:
            md_return = {"{:s}{:s}".format(prefix, key): value for key, value in md_return.items()}

        return md_return

    def align(self, step=0, reflection_angle=0.12, verbosity=3):
        """Align the sample with respect to the beam. GISAXS alignment involves
        vertical translation to the beam center, and rocking theta to get the
        sample plane parralel to the beam. Finally, the angle is re-optimized
        in reflection mode.

        The 'step' argument can optionally be given to jump to a particular
        step in the sequence."""
        start_time = time.time()
        alignment = "Success"
        initial_y = smy.position
        initial_th = sth.position

        align_crazy = self.swing(reflection_angle=reflection_angle)
        crazy_y = smy.position
        crazy_th = sth.position

        if align_crazy[0] == False:
            alignment = "Failed"
            if step <= 4:
                if verbosity >= 4:
                    print("    align: fitting")

                fit_scan(smy, 1.2, 21, fit="HMi")
                ##time.sleep(2)
                fit_scan(sth, 1.5, 21, fit="max")
                # time.sleep(2)

            if step <= 8:
                # fit_scan(smy, 0.3, 21, fit='sigmoid_r')

                fit_edge(smy, 0.6, 21)
                # time.sleep(2)
                # fit_edge(smy, 0.4, 21)
                fit_scan(sth, 0.8, 21, fit="COM")
                # time.sleep(2)
                self.setOrigin(["y", "th"])

            if step <= 9 and reflection_angle is not None:
                # Final alignment using reflected beam
                if verbosity >= 4:
                    print("    align: reflected beam")
                get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2.0)
                # get_beamline().setReflectedBeamROI(total_angle=reflection_angle*2.0, size=[12,2])

                self.thabs(reflection_angle)

                result = fit_scan(sth, 0.4, 41, fit="max")
                # result = fit_scan(sth, 0.2, 81, fit='max') #it's useful for alignment of SmarAct stage
                sth_target = result.values["x_max"] - reflection_angle

                if result.values["y_max"] > 50:
                    th_target = self._axes["th"].motor_to_cur(sth_target)
                    self.thsetOrigin(th_target)

                # fit_scan(smy, 0.2, 21, fit='max')
                # self.setOrigin(['y'])

            if step <= 10:
                self.thabs(0.0)
                beam.off()

        ### save the alignment information
        align_time = time.time() - start_time

        current_data = {
            "a_sample": self.name,
            "i_smx": smx.position,
            "j_smy": smy.position,
            "k_sth": sth.position,
            "l_laserx": laserx.position,
            "m_power": caget(pta.powerV_PV),
            "b_quick_alignment": alignment,
            "c_align_time": align_time,
            "d_offset_y": smy.position - initial_y,
            "e_offset_th": sth.position - initial_th,
            "f_crazy_offset_y": smy.position - crazy_y,
            "g_crazy_offset_th": sth.position - crazy_th,
            "h_search_no": align_crazy[1],
        }

        temp_data = pds.DataFrame([current_data])

        INT_FILENAME = "{}/data/{}.csv".format(os.path.dirname(__file__), "alignment_results.csv")

        if os.path.isfile(INT_FILENAME):
            output_data = pds.read_csv(INT_FILENAME, index_col=0)
            output_data = output_data.append(temp_data, ignore_index=True)
            output_data.to_csv(INT_FILENAME)
        else:
            temp_data.to_csv(INT_FILENAME)

    def align_old(self, step=0, reflection_angle=0.12, verbosity=3):
        """Align the sample with respect to the beam. GISAXS alignment involves
        vertical translation to the beam center, and rocking theta to get the
        sample plane parralel to the beam. Finally, the angle is re-optimized
        in reflection mode.

        The 'step' argument can optionally be given to jump to a particular
        step in the sequence."""
        # if step<=4:
        #     if verbosity>=4:
        #         print('    align: fitting')

        #     fit_scan(smy, 1.2, 21, fit='HMi')
        #     ##time.sleep(2)
        #     fit_scan(sth, 1.5, 21, fit='max')
        ##time.sleep(2)
        cms.modeAlignment()
        self.yo()
        self.tho()
        if step <= 8:
            # fit_scan(smy, 0.3, 21, fit='sigmoid_r')

            fit_edge(smy, 0.6, 21)
            # time.sleep(2)
            # fit_edge(smy, 0.4, 21)
            fit_scan(sth, 0.8, 21, fit="COM")
            # time.sleep(2)
            self.setOrigin(["y", "th"])

        if step <= 9 and reflection_angle is not None:
            # Final alignment using reflected beam
            if verbosity >= 4:
                print("    align: reflected beam")
            get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2.0)
            # get_beamline().setReflectedBeamROI(total_angle=reflection_angle*2.0, size=[12,2])

            self.thabs(reflection_angle)

            result = fit_scan(sth, 0.4, 41, fit="max")
            # result = fit_scan(sth, 0.2, 81, fit='max') #it's useful for alignment of SmarAct stage
            sth_target = result.values["x_max"] - reflection_angle

            if result.values["y_max"] > 50:
                th_target = self._axes["th"].motor_to_cur(sth_target)
                self.thsetOrigin(th_target)

            # fit_scan(smy, 0.2, 21, fit='max')
            # self.setOrigin(['y'])

        if step <= 10:
            self.thabs(0.0)
            beam.off()

    def align_th(self, step=0, reflection_angle=0.12, verbosity=3):
        """Align the sample with respect to the beam. GISAXS alignment involves
        vertical translation to the beam center, and rocking theta to get the
        sample plane parralel to the beam. Finally, the angle is re-optimized
        in reflection mode.

        The 'step' argument can optionally be given to jump to a particular
        step in the sequence."""
        # if step<=4:
        #     if verbosity>=4:
        #         print('    align: fitting')

        #     fit_scan(smy, 1.2, 21, fit='HMi')
        #     ##time.sleep(2)
        #     fit_scan(sth, 1.5, 21, fit='max')
        ##time.sleep(2)
        cms.modeAlignment()
        self.yo()
        self.tho()

        if step <= 9 and reflection_angle is not None:
            # Final alignment using reflected beam
            if verbosity >= 4:
                print("    align: reflected beam")
            get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2.0)
            # get_beamline().setReflectedBeamROI(total_angle=reflection_angle*2.0, size=[12,2])

            self.thabs(reflection_angle)

            result = fit_scan(sth, 0.2, 21, fit="max")
            # result = fit_scan(sth, 0.2, 81, fit='max') #it's useful for alignment of SmarAct stage
            sth_target = result.values["x_max"] - reflection_angle

            if result.values["y_max"] > 50:
                th_target = self._axes["th"].motor_to_cur(sth_target)
                self.thsetOrigin(th_target)

            # fit_scan(smy, 0.2, 21, fit='max')
            # self.setOrigin(['y'])

        if step <= 10:
            self.thabs(0.0)
            beam.off()

    def align_crazy_v2(
        self, step=0, reflection_angle=0.12, ROI_size=[10, 180], th_range=0.3, int_threshold=10, verbosity=3
    ):

        # setting parameters
        rel_th = 1
        ct = 0
        cycle = 0
        intenisty_threshold = 10

        # re-assure the 3 ROI positon
        get_beamline().setDirectBeamROI()
        get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2)

        # set ROI2 as a fixed area
        get_beamline().setROI2ReflectBeamROI(total_angle=reflection_angle * 2, size=ROI_size)
        pilatus2M.roi2.size.y.set(200)
        pilatus2M.roi2.min_xyz.min_y.set(852)

        # def ROI3 in 160pixels with the center located at reflection beam
        # get_beamline().setReflectedBeamROI(total_angle = reflection_angle*2, size=ROI_size) #set ROI3

        # self.thabs(reflection_angle)
        if verbosity >= 4:
            print("  Aligning {}".format(self.name))

        if step <= 0:
            # Prepare for alignment

            if RE.state != "idle":
                RE.abort()

            if get_beamline().current_mode != "alignment":
                # if verbosity>=2:
                # print("WARNING: Beamline is not in alignment mode (mode is '{}')".format(get_beamline().current_mode))
                print("Switching to alignment mode (current mode is '{}')".format(get_beamline().current_mode))
                get_beamline().modeAlignment()

            get_beamline().setDirectBeamROI()

            beam.on()

        if step <= 2:
            if verbosity >= 4:
                print("    align: searching")

            # Estimate full-beam intensity
            value = None
            if True:
                # You can eliminate this, in which case RE.md['beam_intensity_expected'] is used by default
                self.yr(-0.5)
                # detector = gs.DETS[0]
                detector = get_beamline().detector[0]
                value_name = get_beamline().TABLE_COLS[0]
                RE(count([detector]))
                value = detector.read()[value_name]["value"]
                self.yr(0.5)

            if "beam_intensity_expected" in RE.md:
                if value < RE.md["beam_intensity_expected"] * 0.75:
                    print(
                        "WARNING: Direct beam intensity ({}) lower than it should be ({})".format(
                            value, RE.md["beam_intensity_expected"]
                        )
                    )

            # check the last value:
            ii = 0
            while abs(pilatus2M.stats4.total.get() - value) / value < 0.1 and ii < 3:
                ii += 1
                # Find the step-edge
                self.ysearch(
                    step_size=0.2, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1
                )

                # Find the peak
                self.thsearch(step_size=0.2, min_step=0.01, target="max", verbosity=verbosity)

            # last check for height
            self.ysearch(
                step_size=0.05, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1
            )

        # check reflection beam
        self.thr(reflection_angle)
        RE(count([detector]))

        if (
            abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) < 20
            and detector.stats2.max_value.get() > intenisty_threshold
        ):

            # continue the fast alignment
            print("The reflective beam is found! Continue the fast alignment")

            while abs(rel_th) > 0.005 and ct < 5:
                # while detector.roi3.max_value.get() > 50 and ct < 5:

                # absolute beam position
                refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                # roi3 position
                roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                # distance from current postion to the center of roi2 (the disired rel beam position)
                # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                rel_ypos = refl_beam - roi3_beam

                rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                print("The th offset is {}".format(rel_th))
                self.thr(rel_th)

                ct += 1
                RE(count([detector]))

            if detector.stats3.total.get() > 50:

                print("The fast alignment works!")
                self.thr(-reflection_angle)
                self.setOrigin(["y", "th"])

                beam.off()

                return True, ii

            else:
                print("Alignment Error: Cannot Locate the reflection beam")
                self.thr(-reflection_angle)
                beam.off()

                return False, ii

        elif abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) > 5:
            print("Max and Centroid dont Match!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

        else:
            print("Intensiy < threshold!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

    def align_crazy_v3(
        self, step=0, reflection_angle=0.12, ROI_size=[10, 180], th_range=0.3, int_threshold=10, verbosity=3
    ):

        # setting parameters
        rel_th = 1
        ct = 0
        cycle = 0
        intenisty_threshold = 10

        # re-assure the 3 ROI positon
        get_beamline().setDirectBeamROI()
        get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2)
        detector = get_beamline().detector[0]

        # set ROI2 as a fixed area
        get_beamline().setROI2ReflectBeamROI(total_angle=reflection_angle * 2, size=ROI_size)
        pilatus2M.roi2.size.y.set(200)
        pilatus2M.roi2.min_xyz.min_y.set(842)

        # def ROI3 in 160pixels with the center located at reflection beam
        # get_beamline().setReflectedBeamROI(total_angle = reflection_angle*2, size=ROI_size) #set ROI3

        # self.thabs(reflection_angle)
        if verbosity >= 4:
            print("  Aligning {}".format(self.name))

        if step <= 0:
            # Prepare for alignment

            if RE.state != "idle":
                RE.abort()

            if get_beamline().current_mode != "alignment":
                # if verbosity>=2:
                # print("WARNING: Beamline is not in alignment mode (mode is '{}')".format(get_beamline().current_mode))
                print("Switching to alignment mode (current mode is '{}')".format(get_beamline().current_mode))
                get_beamline().modeAlignment()

            get_beamline().setDirectBeamROI()

            beam.on()

        if step <= 2:

            ######################### fast alignment in the case2 and 3 -- NO refl beam
            self.thabs(0.12)
            self.snap(0.5)
            roi2_int = pilatus2M.stats2.total.get()
            roi4_int = pilatus2M.stats4.total.get()
            threshold = 100
            beam_int = 20000
            target_ratio = 0.5
            beam.on()
            if roi2_int < threshold:

                print("CASE 2 or 3")

                roi4_int = pilatus2M.stats4.total.get()
                roi2_int = pilatus2M.stats2.total.get()

                roi4_beam = roi4_int / beam_int

                min_step = 0.005
                # if roi4_beam<target_ratio: #blocking the beam, +Y
                # print(' +Y')
                self.ysearch(
                    step_size=0.01,
                    min_step=0.005,
                    intensity=beam_int,
                    target=0.5,
                    verbosity=verbosity,
                    polarity=-1,
                )
                # else:
                #     print(' -Y')
                #     self.ysearch(step_size=0.01, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)

                roi4_beam = roi4_int / beam_int
                roi2_int = pilatus2M.stats2.total.get()

            # use the beam heigh to find the correct refl beam
            print("Search the refl beam")
            RE(count([pilatus2M]))
            roi4_beam = roi4_int / beam_int
            roi2_int = pilatus2M.stats2.total.get()
            # roi2_int = roi2_i
            th_step = 0.1

            while roi2_int < threshold:
                self.thr(th_step)
                print("th_step {}".format(th_step))
                print("Search the refl beam - th = {}".format(self.thabs()))

                RE(count([pilatus2M]))
                roi4_beam2 = roi4_int / beam_int
                self.yr((roi4_beam2 - roi4_beam) * 0.05)
                if roi4_beam2 < roi4_beam:
                    self.thr(-2 * th_step)
                    self.yr(-2 * (roi4_beam2 - roi4_beam) * 0.05)
                    th_step = -th_step

                    print("REVERSED. th_step {}".format(th_step))

                RE(count([pilatus2M]))
                roi4_beam = roi4_int / beam_int
                roi2_int = pilatus2M.stats2.total.get()

            ######################### fast alignment in the case2 -- y is at 50%

            while abs(rel_th) > 0.005 and ct < 5:
                # while detector.roi3.max_value.get() > 50 and ct < 5:

                print("CASE 2 ")

                # absolute beam position
                refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                # roi3 position
                roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                # distance from current postion to the center of roi2 (the disired rel beam position)
                # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                rel_ypos = refl_beam - roi3_beam

                rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                print("The th offset is {}".format(rel_th))
                self.thr(rel_th)

                ct += 1
                RE(count([pilatus2M]))
                # self.ysearch(step_size=0.01, min_step=0.005, intensity=beam_int, target=0.5, verbosity=verbosity, polarity=-1)

            ######################### fast alignment in the case1 -- both refl and direct beam
            target_ratio = 1
            # self.snap()
            print("CASE 1")

            def get_roi2_4():
                roi2_int = pilatus2M.stats2.total.get()
                roi2_int = roi2_int if roi2_int > 0 else 0
                roi4_int = pilatus2M.stats4.total.get()
                roi4_int = roi4_int if roi4_int > 0 else 0
                return roi2_int / (roi4_int + 10)

            roi2_4 = get_roi2_4()

            min_step = 0.005
            while abs(roi2_4 - target_ratio) > 0.2:
                print(roi2_4)
                if roi2_4 < target_ratio:
                    # print(" +Y")
                    step = min_step
                else:
                    # print(" -Y")
                    step = -min_step

                self.yr(step)
                self.snap()
                roi2_4 = get_roi2_4()

        if step > 5:

            if verbosity >= 4:
                print("    align: searching")

            # Estimate full-beam intensity
            value = None
            if True:
                # You can eliminate this, in which case RE.md['beam_intensity_expected'] is used by default
                self.yr(-0.5)
                # detector = gs.DETS[0]
                detector = get_beamline().detector[0]
                value_name = get_beamline().TABLE_COLS[0]
                RE(count([detector]))
                value = detector.read()[value_name]["value"]
                self.yr(0.5)

            if "beam_intensity_expected" in RE.md:
                if value < RE.md["beam_intensity_expected"] * 0.75:
                    print(
                        "WARNING: Direct beam intensity ({}) lower than it should be ({})".format(
                            value, RE.md["beam_intensity_expected"]
                        )
                    )

            # check the last value:
            ii = 0
            while abs(pilatus2M.stats4.total.get() - value) / value < 0.1 and ii < 3:
                ii += 1
                # Find the step-edge
                self.ysearch(
                    step_size=0.2, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1
                )

                # Find the peak
                self.thsearch(step_size=0.2, min_step=0.01, target="max", verbosity=verbosity)

            # last check for height
            self.ysearch(
                step_size=0.05, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1
            )

        if step > 5:

            # check reflection beam
            self.thr(reflection_angle)
            RE(count([detector]))

            if (
                abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) < 20
                and detector.stats2.max_value.get() > intenisty_threshold
            ):

                # continue the fast alignment
                print("The reflective beam is found! Continue the fast alignment")

                while abs(rel_th) > 0.005 and ct < 5:
                    # while detector.roi3.max_value.get() > 50 and ct < 5:

                    # absolute beam position
                    refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                    # roi3 position
                    roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                    # distance from current postion to the center of roi2 (the disired rel beam position)
                    # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                    rel_ypos = refl_beam - roi3_beam

                    rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                    print("The th offset is {}".format(rel_th))
                    self.thr(rel_th)

                    ct += 1
                    RE(count([detector]))

                # if detector.stats3.total.get()>50:

                #     print('The fast alignment works!')
                #     self.thr(-reflection_angle)
                #     self.setOrigin(['y', 'th'])

                #     beam.off()

                #     return True, ii

                # else:
                #     print('Alignment Error: Cannot Locate the reflection beam')
                #     self.thr(-reflection_angle)
                #     beam.off()

                #     return False, ii

        # elif abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) > 5:
        #     print('Max and Centroid dont Match!')

        #     #perform the full alignment
        #     print('Alignment Error: No reflection beam is found!')
        #     self.thr(-reflection_angle)
        #     beam.off()
        #     return False, ii

        # else:
        #     print('Intensiy < threshold!')

        #     #perform the full alignment
        #     print('Alignment Error: No reflection beam is found!')
        #     self.thr(-reflection_angle)
        #     beam.off()
        #     return False, ii

    def align_crazy_v3_plan(
        self,
        step=0,
        reflection_angle=0.12,
        ROI_size=[10, 180],
        th_range=0.3,
        int_threshold=10,
        direct_beam_int=None,
        verbosity=3,
        detector=None,
        detector_suffix=None,
    ):

        if detector is None:
            # detector = gs.DETS[0]
            detector = get_beamline().detector[0]

        # if detector_suffix is None:
        #     #value_name = gs.TABLE_COLS[0]
        #     value_name = get_beamline().TABLE_COLS[0]
        # else:
        #     value_name = detector.name + detector_suffix

        motors_for_table = [smx, smy, sth]

        @bpp.stage_decorator([detector])
        @bpp.run_decorator(md={})
        @bpp.finalize_decorator(final_plan=shutter_off)
        def inner_align(group=None):
            nonlocal step, reflection_angle

            if group:
                yield from bps.wait(group)

            # setting parameters
            rel_th = 1
            ct = 0
            cycle = 0
            intenisty_threshold = 50

            # re-assure the 3 ROI positon
            get_beamline().setDirectBeamROI()
            get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2)
            detector = get_beamline().detector[0]

            # set ROI2 as a fixed area
            get_beamline().setROI2ReflectBeamROI(total_angle=reflection_angle * 2, size=ROI_size)
            pilatus2M.roi2.size.y.set(200)
            pilatus2M.roi2.min_xyz.min_y.set(842)

            # def ROI3 in 160pixels with the center located at reflection beam
            # get_beamline().setReflectedBeamROI(total_angle = reflection_angle*2, size=ROI_size) #set ROI3

            # self.thabs(reflection_angle)
            if verbosity >= 4:
                print("  Aligning {}".format(self.name))

            if step <= 0:
                print(f"Step <= 0")
                # Prepare for alignment
                if get_beamline().current_mode != "alignment":
                    # if verbosity>=2:
                    # print("WARNING: Beamline is not in alignment mode (mode is '{}')".format(get_beamline().current_mode))
                    print("Switching to alignment mode (current mode is '{}')".format(get_beamline().current_mode))
                    yield from get_beamline().modeAlignment()

                get_beamline().setDirectBeamROI()

                yield from shutter_on()

            if direct_beam_int is not None:
                value = direct_beam_int
            elif hasattr(cms, "direct_beam_int") and cms.direct_beam_int is not None:
                value = cms.direct_beam_int
            else:
                value = 0
                # You can eliminate this, in which case RE.md['beam_intensity_expected'] is used by default
                for n in range(1, 4):
                    self.yr(-0.5)
                    # detector = gs.DETS[0]
                    detector = get_beamline().detector[0]
                    value_name = get_beamline().TABLE_COLS[0]
                    yield from bps.trigger_and_read([detector, *motors_for_table])
                    value = detector.read()[value_name]["value"]
                    if value > 100:
                        cms.direct_beam_int = value
                        self.yr(0.5)
                        break


            if "beam_intensity_expected" in RE.md:
                if value < RE.md["beam_intensity_expected"] * 0.75:
                    print(
                        "WARNING: Direct beam intensity ({}) lower than it should be ({})".format(
                            value, RE.md["beam_intensity_expected"]
                        )
                    )

            if step <= 2:
                print("Step <= 2")

                ######################### fast alignment in the case2 and 3 -- NO refl beam
                self.thabs(0.12)
                # self.snap(0.5)
                yield from bps.trigger_and_read([detector, *motors_for_table])
                roi2_int = pilatus2M.stats2.total.get()
                roi4_int = pilatus2M.stats4.total.get()
                threshold = 500
                beam_int = value
                target_ratio = 0.5
                # yield from shutter_on()
                print(f"roi2_int={roi2_int} threshold={threshold}")

                if roi2_int < threshold:

                    print("CASE 2 or 3")

                    roi4_int = pilatus2M.stats4.total.get()
                    roi2_int = pilatus2M.stats2.total.get()

                    roi4_beam = roi4_int / beam_int

                    min_step = 0.005
                    # if roi4_beam<target_ratio: #blocking the beam, +Y
                    # print(' +Y')
                    # yield from self.search_stub2(
                    #     motor=smy,
                    #     step_size=0.01,
                    #     min_step=0.005,
                    #     target=0.5,
                    #     intensity=beam_int,
                    #     polarity=-1,
                    #     detector=detector,
                    #     detector_suffix="_stats4_total",
                    # )

                    for ii in range(3):
                        norm_stats4 = abs(pilatus2M.stats4.total.get() - beam_int) / beam_int
                        print(f"ii={ii} norm_stats4={norm_stats4}")

                        if ii > 0 and norm_stats4 > 0.2: 
                            break

                        # Find the step-edge
                        yield from self.search_stub2(
                            motor=smy,
                            step_size=0.2,
                            min_step=0.005,
                            target=0.5,
                            intensity=beam_int,
                            polarity=-1,
                            detector=detector,
                            detector_suffix="_stats4_total",
                        )

                        yield from self.search_stub2(
                            motor=sth,
                            step_size=0.2,
                            min_step=0.01,
                            target="max",
                            polarity=-1,
                            detector=detector,
                            detector_suffix="_stats4_total",
                        )

                        # self.ysearch(step_size=0.2, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)

                        # # Find the peak
                        # self.thsearch(step_size=0.2, min_step=0.01, target='max', verbosity=verbosity)

                    # last check for height
                    # self.ysearch(step_size=0.05, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)
                    yield from self.search_stub2(
                        motor=smy,
                        step_size=0.05,
                        min_step=0.005,
                        target=0.5,
                        intensity=beam_int,
                        polarity=-1,
                        detector=detector,
                        detector_suffix="_stats4_total",
                    )

                    # self.ysearch(step_size=0.01, min_step=0.005, intensity=beam_int, target=0.5, verbosity=verbosity, polarity=-1)
                    # else:
                    #     print(' -Y')
                    #     self.ysearch(step_size=0.01, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)

                    roi4_beam = roi4_int / beam_int
                    roi2_int = pilatus2M.stats2.total.get()

                else:
                    #very close to aligned position
                    reflection_angle = 0
                # #use the beam heigh to find the correct refl beam
                # print('Search the refl beam')
                # yield from bps.trigger_and_read([detector, *motors_for_table])
                # # RE(count([pilatus2M]))

                # roi4_beam = roi4_int/beam_int
                # roi2_int = pilatus2M.stats2.total.get()
                # # roi2_int = roi2_i
                # th_step = 0.1

                # while roi2_int<threshold:
                #     self.thr(th_step)
                #     print('th_step {}'.format(th_step))
                #     print('Search the refl beam - th = {}'.format(self.thabs()))

                #     yield from bps.trigger_and_read([detector, *motors_for_table])
                #     # RE(count([pilatus2M]))
                #     roi4_beam2 = roi4_int/beam_int
                #     self.yr((roi4_beam2-roi4_beam)*0.05)
                #     if roi4_beam2 < roi4_beam:
                #         self.thr(-2*th_step)
                #         self.yr(-2*(roi4_beam2-roi4_beam)*0.05)
                #         th_step = -th_step

                #         print('REVERSED. th_step {}'.format(th_step))

                #     yield from bps.trigger_and_read([detector, *motors_for_table])
                #     # RE(count([pilatus2M]))
                #     roi4_beam = roi4_int/beam_int
                #     roi2_int = pilatus2M.stats2.total.get()

                # ######################### fast alignment in the case2 -- y is at 50%

                # while abs(rel_th) > 0.005 and ct < 5:
                # # while detector.roi3.max_value.get() > 50 and ct < 5:

                #     print("CASE 2 ")

                #     #absolute beam position
                #     refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                #     #roi3 position
                #     roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get()/2

                #     #distance from current postion to the center of roi2 (the disired rel beam position)
                #     # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                #     rel_ypos = refl_beam - roi3_beam

                #     rel_th = rel_ypos/get_beamline().SAXS.distance/1000*0.172/np.pi*180/2

                #     print('The th offset is {}'.format(rel_th))
                #     self.thr(rel_th)

                #     ct += 1
                #     yield from bps.trigger_and_read([detector, *motors_for_table])
                #     # RE(count([pilatus2M]))
                #     # self.ysearch(step_size=0.01, min_step=0.005, intensity=beam_int, target=0.5, verbosity=verbosity, polarity=-1)

                ######################### fast alignment in the case1 -- both refl and direct beam

                # self.thr(reflection_angle)
                # yield from bps.trigger_and_read([detector, *motors_for_table])
                # RE(count([detector]))

            #     # if abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) < 20 and detector.stats2.max_value.get() > intenisty_threshold:

            #     target_ratio = 1
            #     # self.snap()
            #     print("CASE 1")

            #     def get_roi2_4():
            #         roi2_int = pilatus2M.stats2.total.get()
            #         roi2_int = roi2_int if roi2_int > 0 else 0
            #         roi4_int = pilatus2M.stats4.total.get()
            #         roi4_int = roi4_int if roi4_int > 0 else 0
            #         return roi2_int/(roi4_int + 10)

            #     roi2_4 = get_roi2_4()

            #     min_step=0.005
            #     while abs(roi2_4 - target_ratio)>0.2:
            #         print(roi2_4)
            #         if roi2_4<target_ratio:
            #             print(' +Y')
            #             step = min_step
            #         else:
            #             print(' -Y')
            #             step = -min_step

            #         self.yr(step)
            #         yield from bps.trigger_and_read([detector, *motors_for_table])
            #         # self.snap()
            #         roi2_4 = get_roi2_4()

            # if step>5:

            #     if verbosity>=4:
            #         print('    align: searching')

            #     # Estimate full-beam intensity
            #     value = None
            #     if True:
            #         # You can eliminate this, in which case RE.md['beam_intensity_expected'] is used by default
            #         self.yr(-0.5)
            #         #detector = gs.DETS[0]
            #         detector = get_beamline().detector[0]
            #         value_name = get_beamline().TABLE_COLS[0]
            #         yield from bps.trigger_and_read([detector, *motors_for_table])
            #         # RE(count([detector]))
            #         value = detector.read()[value_name]['value']
            #         self.yr(0.5)

            #     if 'beam_intensity_expected' in RE.md:
            #         if value<RE.md['beam_intensity_expected']*0.75:
            #             print('WARNING: Direct beam intensity ({}) lower than it should be ({})'.format(value, RE.md['beam_intensity_expected']))

            #     #check the last value:
            #     ii = 0
            #     while abs(pilatus2M.stats4.total.get() - value)/value < 0.1 and ii < 3:
            #         ii += 1
            #         # Find the step-edge
            #         yield from self.search_stub2(
            #             motor=smy, step_size=0.2, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity,
            #             polarity=-1, detector=detector, detector_suffix='_stats4_total'
            #         )
            #         # self.ysearch(step_size=0.2, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)

            #         # Find the peak
            #         yield from self.search_stub2(
            #             motor=sth, step_size=0.2, min_step=0.01, target='max', verbosity=verbosity, detector=detector,
            #             detector_suffix='_stats4_total'
            #         )
            #         # self.thsearch(step_size=0.2, min_step=0.01, target='max', verbosity=verbosity)

            #     #last check for height
            #     yield from self.search_stub2(
            #         motor=smy, step_size=0.05, min_step=0.005, intensity=value, target=0.5,
            #         verbosity=verbosity, polarity=-1, detector=detector, detector_suffix='_stats4_total'
            #     )
            #     # self.ysearch(step_size=0.05, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)

            if step < 5:
                print(f"Step <= 5")

                # check reflection beam
                self.thr(reflection_angle)
                yield from bps.trigger_and_read([detector, *motors_for_table])
                # RE(count([detector]))

                stat2_max_xy_centr = abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y)
                stat2_max_value = detector.stats2.max_value.get()
                print(f"stat2_max_xy_centr={stat2_max_xy_centr} stat2_max_value={stat2_max_value}")

                if (stat2_max_xy_centr < 20 and stat2_max_value > intenisty_threshold):

                    # continue the fast alignment
                    print("The reflective beam is found! Continue the fast alignment")

                    while abs(rel_th) > 0.005 and ct < 5:
                        # while detector.roi3.max_value.get() > 50 and ct < 5:

                        # absolute beam position
                        refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                        # roi3 position
                        roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                        # distance from current postion to the center of roi2 (the disired rel beam position)
                        # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                        rel_ypos = refl_beam - roi3_beam

                        rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                        print("The th offset is {}".format(rel_th))
                        self.thr(rel_th)

                        ct += 1
                        yield from bps.trigger_and_read([detector, *motors_for_table])
                        # RE(count([detector]))

                    # if detector.stats3.total.get()>50:

                    #     print('The fast alignment works!')
                    #     self.thr(-reflection_angle)
                    #     self.setOrigin(['y', 'th'])

                    #     beam.off()

                    #     return True, ii

                    # else:
                    #     print('Alignment Error: Cannot Locate the reflection beam')
                    #     self.thr(-reflection_angle)
                    #     beam.off()

                    #     return False, ii

            # elif abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) > 5:
            #     print('Max and Centroid dont Match!')

            #     #perform the full alignment
            #     print('Alignment Error: No reflection beam is found!')
            #     self.thr(-reflection_angle)
            #     beam.off()
            #     return False, ii

            # else:
            #     print('Intensiy < threshold!')

            #     #perform the full alignment
            #     print('Alignment Error: No reflection beam is found!')
            #     self.thr(-reflection_angle)
            #     beam.off()
            #     return False, ii

        group_name = "setup_aligment"

        #alignment mode
        yield from bps.abs_set(bsx, cms.bsx_pos + 3, group=group_name)
        beam.setTransmission(1e-6)

        #align as abovve
        yield from inner_align(group=group_name)

        #move bs back
        yield from bps.abs_set(bsx, cms.bsx_pos, group=group_name)
        yield from bps.wait(group_name)

        #set the position for sample
        # self.thr(reflection_angle)
        # self.setOrigin(['y', 'th'])

    def swing_v2(
        self, step=0, reflection_angle=0.12, ROI_size=[10, 180], th_range=0.3, int_threshold=10, verbosity=3
    ):

        # setting parameters
        rel_th = 1
        ct = 0
        cycle = 0
        intenisty_threshold = 10

        # re-assure the 3 ROI positon
        get_beamline().setDirectBeamROI()
        get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2)

        # set ROI2 as a fixed area
        get_beamline().setROI2ReflectBeamROI(total_angle=reflection_angle * 2, size=ROI_size)
        pilatus2M.roi2.size.y.set(190)
        pilatus2M.roi2.min_xyz.min_y.set(852)

        # def ROI3 in 160pixels with the center located at reflection beam
        # get_beamline().setReflectedBeamROI(total_angle = reflection_angle*2, size=ROI_size) #set ROI3

        # self.thabs(reflection_angle)
        if verbosity >= 4:
            print("  Aligning {}".format(self.name))

            # if step<=0:
            #     # Prepare for alignment

            #     if RE.state!='idle':
            #         RE.abort()

            #     if get_beamline().current_mode!='alignment':
            #         #if verbosity>=2:
            #             #print("WARNING: Beamline is not in alignment mode (mode is '{}')".format(get_beamline().current_mode))
            #         print("Switching to alignment mode (current mode is '{}')".format(get_beamline().current_mode))
            #         get_beamline().modeAlignment()

            get_beamline().setDirectBeamROI()

            beam.on()

        if step <= 2:
            # if verbosity>=4:
            #     print('    align: searching')

            # Estimate full-beam intensity
            value = None
            if True:
                # You can eliminate this, in which case RE.md['beam_intensity_expected'] is used by default
                self.yr(-0.5)
                # detector = gs.DETS[0]
                detector = get_beamline().detector[0]
                # value_name = get_beamline().TABLE_COLS[0]
                beam.on()
                RE(count([detector]))
                value = detector.read()["pilatus2M_stats4_total"]["value"]
                self.yr(0.5)

            # if 'beam_intensity_expected' in RE.md:
            #     if value<RE.md['beam_intensity_expected']*0.75:
            #         print('WARNING: Direct beam intensity ({}) lower than it should be ({})'.format(value, RE.md['beam_intensity_expected']))

            # check the last value:
            # value=20000
            ii = 0
            while abs(pilatus2M.stats4.total.get() - value) / value < 0.1 and ii < 3:
                ii += 1
                # Find the step-edge
                fastsearch = RE(
                    self.search_plan(
                        motor=smy,
                        step_size=0.1,
                        min_step=0.01,
                        target=0.5,
                        intensity=20000,
                        polarity=-1,
                        fastsearch=True,
                        detector_suffix="_stats4_total",
                    )
                )
                if fastsearch == True:
                    break
                # Find the peak
                # self.thsearch(step_size=0.2, min_step=0.01, target='max', verbosity=verbosity)
                fastsearch = RE(
                    self.search_plan(
                        motor=sth,
                        step_size=0.2,
                        min_step=0.01,
                        target="max",
                        fastsearch=True,
                        detector_suffix="_stats4_total",
                    )
                )
                if fastsearch == True:
                    break
            # last check for height
            # self.ysearch(step_size=0.05, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)
            if fastsearch == False:
                RE(
                    self.search_plan(
                        motor=smy,
                        step_size=0.05,
                        min_step=0.005,
                        target=0.5,
                        intensity=20000,
                        polarity=-1,
                        detector_suffix="_stats4_total",
                    )
                )

        # check reflection beam
        self.thr(reflection_angle)
        RE(count([detector]))

        if (
            abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) < 20
            and detector.stats2.max_value.get() > intenisty_threshold
        ):

            # continue the fast alignment
            print("The reflective beam is found! Continue the fast alignment")

            while abs(rel_th) > 0.005 and ct < 5:
                # while detector.roi3.max_value.get() > 50 and ct < 5:

                # absolute beam position
                refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                # roi3 position
                roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                # distance from current postion to the center of roi2 (the disired rel beam position)
                # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                rel_ypos = refl_beam - roi3_beam

                rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                print("The th offset is {}".format(rel_th))
                self.thr(rel_th)

                ct += 1
                RE(count([detector]))

            if detector.stats3.total.get() > 50:

                print("The fast alignment works!")
                self.thr(-reflection_angle)

                if fastsearch == False:
                    RE(
                        self.search_plan(
                            motor=smy,
                            step_size=0.05,
                            min_step=0.005,
                            target=0.5,
                            intensity=20000,
                            polarity=-1,
                            detector_suffix="_stats4_total",
                        )
                    )

                self.setOrigin(["y", "th"])

                beam.off()

                return True, ii

            else:
                print("Alignment Error: Cannot Locate the reflection beam")
                self.thr(-reflection_angle)
                beam.off()

                return False, ii

        elif abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) > 5:
            print("Max and Centroid dont Match!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

        else:
            print("Intensiy < threshold!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

    def swing_March(
        self, step=0, reflection_angle=0.12, ROI_size=[10, 180], th_range=0.3, intensity=20000, int_threshold=10, verbosity=3
    ):

        # setting parameters
        rel_th = 1
        ct = 0
        cycle = 0
        intenisty_threshold = 10

        # re-assure the 3 ROI positon
        get_beamline().setDirectBeamROI()
        get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2)

        # set ROI2 as a fixed area
        get_beamline().setROI2ReflectBeamROI(total_angle=reflection_angle * 2, size=ROI_size)
        pilatus2M.roi2.size.y.set(190)
        pilatus2M.roi2.min_xyz.min_y.set(852)

        # def ROI3 in 160pixels with the center located at reflection beam
        # get_beamline().setReflectedBeamROI(total_angle = reflection_angle*2, size=ROI_size) #set ROI3

        # self.thabs(reflection_angle)
        if verbosity >= 4:
            print("  Aligning {}".format(self.name))

            # if step<=0:
            #     # Prepare for alignment

            #     if RE.state!='idle':
            #         RE.abort()

            #     if get_beamline().current_mode!='alignment':
            #         #if verbosity>=2:
            #             #print("WARNING: Beamline is not in alignment mode (mode is '{}')".format(get_beamline().current_mode))
            #         print("Switching to alignment mode (current mode is '{}')".format(get_beamline().current_mode))
            #         get_beamline().modeAlignment()

            get_beamline().setDirectBeamROI()

            beam.on()

        if step <= 2:
            # if verbosity>=4:
            #     print('    align: searching')

            # Estimate full-beam intensity
            value = None
            if True:
                # You can eliminate this, in which case RE.md['beam_intensity_expected'] is used by default
                # self.yr(-0.5)
                # detector = gs.DETS[0]
                detector = get_beamline().detector[0]
                # value_name = get_beamline().TABLE_COLS[0]
                beam.on()
                RE(count([detector]))
                value = detector.read()["pilatus2M_stats4_total"]["value"]
                # self.yr(0.5)

            # if 'beam_intensity_expected' in RE.md:
            #     if value<RE.md['beam_intensity_expected']*0.75:
            #         print('WARNING: Direct beam intensity ({}) lower than it should be ({})'.format(value, RE.md['beam_intensity_expected']))

            # check the last value:
            value=20000
            ii = 0
            while abs(pilatus2M.stats4.total.get() - value) / value < 0.1 and ii < 3:
                ii += 1
                # Find the step-edge
                fastsearch = RE(
                    self.search_plan(
                        motor=smy,
                        step_size=0.1,
                        min_step=0.01,
                        target=0.5,
                        intensity=20000,
                        polarity=-1,
                        # fastsearch=True,
                        detector_suffix="_stats4_total",
                    )
                )
                if fastsearch == True:
                    break
            #     # Find the peak
            #     # self.thsearch(step_size=0.2, min_step=0.01, target='max', verbosity=verbosity)
            #     fastsearch = RE(
            #         self.search_plan(
            #             motor=sth,
            #             step_size=0.2,
            #             min_step=0.01,
            #             target="max",
            #             fastsearch=True,
            #             detector_suffix="_stats4_total",
            #         )
            #     )
            #     if fastsearch == True:
            #         break
            # # last check for height
            # # self.ysearch(step_size=0.05, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)
            # if fastsearch == False:
            #     RE(
            #         self.search_plan(
            #             motor=smy,
            #             step_size=0.05,
            #             min_step=0.005,
            #             target=0.5,
            #             intensity=20000,
            #             polarity=-1,
            #             detector_suffix="_stats4_total",
            #         )
            #     )

        # check reflection beam
        self.thr(reflection_angle)
        RE(count([detector]))

        if (
            abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) < 20
            and detector.stats2.max_value.get() > intenisty_threshold
        ):

            # continue the fast alignment
            print("The reflective beam is found! Continue the fast alignment")

            #for sth
            while abs(rel_th) > 0.005 and ct < 5:
                # while detector.roi3.max_value.get() > 50 and ct < 5:

                # absolute beam position
                refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                # roi3 position
                roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                # distance from current postion to the center of roi2 (the disired rel beam position)
                # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                rel_ypos = refl_beam - roi3_beam

                rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                print("The th offset is {}".format(rel_th))
                self.thr(rel_th)

                ct += 1
                RE(count([detector]))

            #for smy
            # Find the step-edge
            fastsearch = RE(
                self.search_plan(
                    motor=smy,
                    step_size=0.05,
                    min_step=0.01,
                    target="max",
                    # intensity=intensity,
                    polarity=-1,
                    # fastsearch=True,
                    detector_suffix="_stats3_total",
                )
            )


            # if detector.stats3.total.get() > 50:

            #     print("The fast alignment works!")
            #     self.thr(-reflection_angle)

            #     if fastsearch == False:
            #         RE(
            #             self.search_plan(
            #                 motor=smy,
            #                 step_size=0.05,
            #                 min_step=0.005,
            #                 target=0.5,
            #                 intensity=20000,
            #                 polarity=-1,
            #                 detector_suffix="_stats4_total",
            #             )
            #         )

            self.setOrigin(["y", "th"])

            beam.off()

            return True, ii

            # else:
            #     print("Alignment Error: Cannot Locate the reflection beam")
            #     self.thr(-reflection_angle)
            #     beam.off()

            #     return False, ii

        # elif abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) > 5:
        #     print("Max and Centroid dont Match!")

        #     # perform the full alignment
        #     print("Alignment Error: No reflection beam is found!")
        #     self.thr(-reflection_angle)
        #     beam.off()
        #     return False, ii

        # else:
        #     print("Intensiy < threshold!")

        #     # perform the full alignment
        #     print("Alignment Error: No reflection beam is found!")
        #     self.thr(-reflection_angle)
        #     beam.off()
        #     return False, ii
        
    def swing(
        self, step=0, reflection_angle=0.12, ROI_size=[10, 180], th_range=0.3, int_threshold=10, verbosity=3
    ):

        # setting parameters
        rel_th = 1
        ct = 0
        cycle = 0
        intenisty_threshold = 10

        # re-assure the 3 ROI positon
        get_beamline().setDirectBeamROI()
        get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2)

        # set ROI2 as a fixed area
        get_beamline().setROI2ReflectBeamROI(total_angle=reflection_angle * 2, size=ROI_size)
        pilatus2M.roi2.size.y.set(190)
        pilatus2M.roi2.min_xyz.min_y.set(852)

        # def ROI3 in 160pixels with the center located at reflection beam
        # get_beamline().setReflectedBeamROI(total_angle = reflection_angle*2, size=ROI_size) #set ROI3

        # self.thabs(reflection_angle)
        if verbosity >= 4:
            print("  Aligning {}".format(self.name))

            # if step<=0:
            #     # Prepare for alignment

            #     if RE.state!='idle':
            #         RE.abort()

            #     if get_beamline().current_mode!='alignment':
            #         #if verbosity>=2:
            #             #print("WARNING: Beamline is not in alignment mode (mode is '{}')".format(get_beamline().current_mode))
            #         print("Switching to alignment mode (current mode is '{}')".format(get_beamline().current_mode))
            #         get_beamline().modeAlignment()

            get_beamline().setDirectBeamROI()

            RE(beam.on())

        if step <= 2:
            # if verbosity>=4:
            #     print('    align: searching')

            # Estimate full-beam intensity
            value = None
            if True:
                # You can eliminate this, in which case RE.md['beam_intensity_expected'] is used by default
                self.yr(-1)
                # detector = gs.DETS[0]
                detector = get_beamline().detector[0]
                # value_name = get_beamline().TABLE_COLS[0]
                RE(beam.on())
                RE(count([detector]))
                value = detector.read()["pilatus2M_stats4_total"]["value"]
                self.yr(1)

            # if 'beam_intensity_expected' in RE.md:
            #     if value<RE.md['beam_intensity_expected']*0.75:
            #         print('WARNING: Direct beam intensity ({}) lower than it should be ({})'.format(value, RE.md['beam_intensity_expected']))

            # check the last value:
            # value=20000
            ii = 0
            while abs(pilatus2M.stats4.total.get() - value) / value < 0.1 and ii < 3:
                ii += 1
                # Find the step-edge
                RE(
                    self.search_plan(
                        motor=smy,
                        step_size=0.1,
                        min_step=0.01,
                        target=0.5,
                        intensity=20000,
                        polarity=-1,
                        detector_suffix="_stats4_total",
                    )
                )
                # Find the peak
                # self.thsearch(step_size=0.2, min_step=0.01, target='max', verbosity=verbosity)
                RE(
                    self.search_plan(
                        motor=sth, step_size=0.2, min_step=0.01, target="max", detector_suffix="_stats4_total"
                    )
                )
            # last check for height
            # self.ysearch(step_size=0.05, min_step=0.005, intensity=value, target=0.5, verbosity=verbosity, polarity=-1)
            RE(
                self.search_plan(
                    motor=smy,
                    step_size=0.05,
                    min_step=0.005,
                    target=0.5,
                    intensity=20000,
                    polarity=-1,
                    detector_suffix="_stats4_total",
                )
            )

        # check reflection beam
        self.thr(reflection_angle)
        RE(count([detector]))

        if (
            abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) < 20
            and detector.stats2.max_value.get() > intenisty_threshold
        ):

            # continue the fast alignment
            print("The reflective beam is found! Continue the fast alignment")

            while abs(rel_th) > 0.005 and ct < 5:
                # while detector.roi3.max_value.get() > 50 and ct < 5:

                # absolute beam position
                refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                # roi3 position
                roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                # distance from current postion to the center of roi2 (the disired rel beam position)
                # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                rel_ypos = refl_beam - roi3_beam

                rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                print("The th offset is {}".format(rel_th))
                self.thr(rel_th)

                ct += 1
                RE(count([detector]))

            if detector.stats3.total.get() > 50:

                print("The fast alignment works!")
                self.thr(-reflection_angle)

                self.setOrigin(["y", "th"])

                beam.off()

                return True, ii

            else:
                print("Alignment Error: Cannot Locate the reflection beam")
                self.thr(-reflection_angle)
                beam.off()

                return False, ii

        elif abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) > 5:
            print("Max and Centroid dont Match!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

        else:
            print("Intensiy < threshold!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

    def crazy_th(
        self, step=0, reflection_angle=0.12, ROI_size=[10, 180], th_range=0.3, int_threshold=10, verbosity=3
    ):

        # setting parameters
        rel_th = 1
        ct = 0
        cycle = 0
        intenisty_threshold = 10

        # re-assure the 3 ROI positon
        get_beamline().setDirectBeamROI()
        get_beamline().setReflectedBeamROI(total_angle=reflection_angle * 2)

        # set ROI2 as a fixed area
        get_beamline().setROI2ReflectBeamROI(total_angle=reflection_angle * 2, size=ROI_size)
        pilatus2M.roi2.size.y.set(200)
        pilatus2M.roi2.min_xyz.min_y.set(852)

        # def ROI3 in 160pixels with the center located at reflection beam
        # get_beamline().setReflectedBeamROI(total_angle = reflection_angle*2, size=ROI_size) #set ROI3

        # self.thabs(reflection_angle)
        if verbosity >= 4:
            print("  Aligning {}".format(self.name))

        if step <= 0:
            # Prepare for alignment

            if RE.state != "idle":
                RE.abort()

            if get_beamline().current_mode != "alignment":
                # if verbosity>=2:
                # print("WARNING: Beamline is not in alignment mode (mode is '{}')".format(get_beamline().current_mode))
                print("Switching to alignment mode (current mode is '{}')".format(get_beamline().current_mode))
                get_beamline().modeAlignment()

            get_beamline().setDirectBeamROI()

            beam.on()

        detector = pilatus2M
        RE(pilatus2M.setExposureTime(0.5))
        self.thabs(reflection_angle)
        RE(count([detector]))

        if (
            abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) < 20
            and detector.stats2.max_value.get() > intenisty_threshold
        ):

            # continue the fast alignment
            print("The reflective beam is found! Continue the fast alignment")

            while abs(rel_th) > 0.005 and ct < 5:
                # while detector.roi3.max_value.get() > 50 and ct < 5:

                # absolute beam position
                refl_beam = detector.roi2.min_xyz.min_y.get() + detector.stats2.max_xy.y.get()

                # roi3 position
                roi3_beam = detector.roi3.min_xyz.min_y.get() + detector.roi3.size.y.get() / 2

                # distance from current postion to the center of roi2 (the disired rel beam position)
                # rel_ypos = detector.stats2.max_xy.get().y - detector.stats2.size.get().y
                rel_ypos = refl_beam - roi3_beam

                rel_th = rel_ypos / get_beamline().SAXS.distance / 1000 * 0.172 / np.pi * 180 / 2

                print("The th offset is {}".format(rel_th))
                self.thr(rel_th)

                ct += 1
                RE(count([detector]))

            if detector.stats3.total.get() > 50:

                print("The fast alignment works!")
                self.thr(-reflection_angle)
                self.setOrigin(["y", "th"])

                beam.off()

                return True, ii

            else:
                print("Alignment Error: Cannot Locate the reflection beam")
                self.thr(-reflection_angle)
                beam.off()

                return False, ii

        elif abs(detector.stats2.max_xy.get().y - detector.stats2.centroid.get().y) > 5:
            print("Max and Centroid dont Match!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

        else:
            print("Intensiy < threshold!")

            # perform the full alignment
            print("Alignment Error: No reflection beam is found!")
            self.thr(-reflection_angle)
            beam.off()
            return False, ii

    def measureInitial(self, exposure_time=10, bounds=[0, 50]):

        pos_list = np.meshgrid(bounds[0], bounds[1], 2)
        # pos_list.append([np.average(bounds),np.average(bounds)])
        command["out_of_bound"] = False

        for x_pos, y_pos in pos_list:
            start_time = time.time()
            if verbosity >= 3:
                print(
                    "{}Driving to point {}/{}; (x,yy) = ({:.3f}, {:.3f})".format(
                        prefix, imeasure, num_to_measure, x_pos, yy_pos
                    )
                )

            self.xabs(x_pos)
            # self.yabs(y_pos)
            self.yyabs(yy_pos)

            while smx.moving == True or smy2.moving == True:
                time.sleep(1)
            while abs(self.xpos(verbosity=0) - x_pos) > 0.1 or abs(self.yypos(verbosity=0) - yy_pos) > 0.1:
                time.sleep(1)

            self.measure(exposure_time=exposure_time, extra=extra, **md)
            header = db[-1]  # The most recent measurement
            # command['filename'] = '{}'.format(header.start['filename'][:-1])
            command["filename"] = "{}".format(header.start["filename"])
            command["x_position"] = self.xpos(verbosity=0)
            command["y_position"] = self.yypos(verbosity=0)
            command["h1_position"] = self.xy2h(self.xpos(), self.yypos())[0]
            command["h2_position"] = self.xy2h(self.xpos(), self.yypos())[1]

            cost_time = time.time() - start_time

            command["cost"] = cost_time

            command["h1_para"] = self.para1
            command["h2_para"] = self.para2

            command["measured"] = True
            command["analyzed"] = False

        measure_queue.publish(commands)  # Send results for analysis

    def search_plan(
        self,
        motor=smy,
        step_size=0.2,
        min_step=0.05,
        intensity=None,
        target=0.5,
        detector=None,
        detector_suffix=None,
        polarity=-1,
        verbosity=3,
    ):
        """Moves this axis, searching for a target value.

        Parameters
        ----------
        step_size : float
            The initial step size when moving the axis
        min_step : float
            The final (minimum) step size to try
        intensity : float
            The expected full-beam intensity readout
        target : 0.0 to 1.0
            The target ratio of full-beam intensity; 0.5 searches for half-max.
            The target can also be 'max' to find a local maximum.
        detector, detector_suffix
            The beamline detector (and suffix, such as '_stats4_total') to trigger to measure intensity
        polarity : +1 or -1
            Positive motion assumes, e.g. a step-height 'up' (as the axis goes more positive)
        """
        print("HERE!!")

        if detector is None:
            # detector = gs.DETS[0]
            detector = get_beamline().detector[0]
        if detector_suffix is None:
            # value_name = gs.TABLE_COLS[0]
            value_name = get_beamline().TABLE_COLS[0]
        else:
            value_name = detector.name + detector_suffix

        print(f"detector={detector}")

        @bpp.stage_decorator([detector])
        @bpp.run_decorator(md={})
        def inner_search():
            nonlocal intensity, target, step_size

            if not get_beamline().beam.is_on():
                print("WARNING: Experimental shutter is not open.")

            if intensity is None:
                intensity = RE.md["beam_intensity_expected"]

            # bec.disable_table()

            # Check current value
            vv = yield from bps.trigger_and_read([detector, motor])
            value = vv[value_name]["value"]
            # RE(count([detector]))
            # value = detector.read()[value_name]['value']

            if target == "max":

                if verbosity >= 5:
                    print("Performing search on axis '{}' target is 'max'".format(self.name))

                max_value = value
                # max_position = self.get_position(verbosity=0)

                direction = +1 * polarity

                while step_size >= min_step:
                    if verbosity >= 4:
                        print("        move {} by {} × {}".format(self.name, direction, step_size))

                    #  pos = yield from bps.rd(motor)
                    yield from bps.mvr(motor, direction * step_size)
                    # self.move_relative(move_amount=direction*step_size, verbosity=verbosity-2)

                    prev_value = value
                    yield from bps.trigger_and_read([detector, motor])
                    # RE(count([detector]))

                    value = detector.read()[value_name]["value"]
                    # if verbosity>=3:
                    #     print("      {} = {:.3f} {}; value : {}".format(self.name, self.get_position(verbosity=0), self.units, value))

                    if value > max_value:
                        max_value = value
                        # max_position = self.get_position(verbosity=0)

                    if value > prev_value:
                        # Keep going in this direction...
                        pass
                    else:
                        # Switch directions!
                        direction *= -1
                        step_size *= 0.5

            elif target == "min":

                if verbosity >= 5:
                    print("Performing search on axis '{}' target is 'min'".format(self.name))

                direction = +1 * polarity

                while step_size >= min_step:
                    if verbosity >= 4:
                        print("        move {} by {} × {}".format(self.name, direction, step_size))

                    # pos = yield from bps.rd(motor)
                    yield from bps.mvr(motor, direction * step_size)
                    # self.move_relative(move_amount=direction*step_size, verbosity=verbosity-2)

                    prev_value = value
                    yield from bps.trigger_and_read([detector, motor])
                    # RE(count([detector]))
                    value = detector.read()[value_name]["value"]
                    if verbosity >= 3:
                        print(
                            "      {} = {:.3f} {}; value : {}".format(
                                self.name, self.get_position(verbosity=0), self.units, value
                            )
                        )

                    if value < prev_value:
                        # Keep going in this direction...
                        pass
                    else:
                        # Switch directions!
                        direction *= -1
                        step_size *= 0.5

            else:

                target_rel = target
                target = target_rel * intensity

                if verbosity >= 5:
                    print(
                        "Performing search on axis '{}' target {} × {} = {}".format(
                            self.name, target_rel, intensity, target
                        )
                    )
                if verbosity >= 4:
                    print("      value : {} ({:.1f}%)".format(value, 100.0 * value / intensity))

                # Determine initial motion direction
                if value > target:
                    direction = -1 * polarity
                else:
                    direction = +1 * polarity

                while step_size >= min_step:

                    if verbosity >= 4:
                        print("        move {} by {} × {}".format(self.name, direction, step_size))

                    # pos = yield from bps.rd(motor)
                    yield from bps.mvr(motor, direction * step_size)
                    # self.move_relative(move_amount=direction*step_size, verbosity=verbosity-2)

                    yield from bps.trigger_and_read([detector, motor])
                    # RE(count([detector]))
                    value = detector.read()[value_name]["value"]
                    # if verbosity>=3:
                    #    print("      {} = {:.3f} {}; value : {} ({:.1f}%)".format(self.name, self.get_position(verbosity=0), self.units, value, 100.0*value/intensity))

                    # Determine direction
                    if value > target:
                        new_direction = -1.0 * polarity
                    else:
                        new_direction = +1.0 * polarity

                    if abs(direction - new_direction) < 1e-4:
                        # Same direction as we've been going...
                        # ...keep moving this way
                        pass
                    else:
                        # Switch directions!
                        direction *= -1
                        step_size *= 0.5

            # bec.enable_table()

        yield from inner_search()

    def search_stub2(
        self,
        motor=smy,
        step_size=1.0,
        min_step=0.05,
        intensity=None,
        target=0.5,
        detector=None,
        detector_suffix=None,
        polarity=1,
        verbosity=3,
    ):

        if detector is None:
            # detector = gs.DETS[0]
            detector = get_beamline().detector[0]
        if detector_suffix is None:
            # value_name = gs.TABLE_COLS[0]
            value_name = get_beamline().TABLE_COLS[0]
        else:
            value_name = detector.name + detector_suffix

        if intensity is None:
            intensity = RE.md["beam_intensity_expected"]

        motors_for_table = [smx, smy, sth]

        # Check current value
        vv = yield from bps.trigger_and_read([detector, *motors_for_table])
        value = vv[value_name]["value"]

        if target == "max":

            if verbosity >= 5:
                print("Performing search on axis '{}' target is 'max'".format(self.name))

            max_value = value
            # max_position = self.get_position(verbosity=0)

            direction = polarity

            while step_size >= min_step:
                if verbosity >= 4:
                    print("        move {} by {} × {}".format(self.name, direction, step_size))

                yield from bps.mvr(motor, direction * step_size)

                prev_value = value
                yield from bps.trigger_and_read([detector, *motors_for_table])

                value = detector.read()[value_name]["value"]
                # if verbosity>=3:
                #     print("      {} = {:.3f} {}; value : {}".format(self.name, self.get_position(verbosity=0), self.units, value))

                if value > max_value:
                    max_value = value
                    # max_position = self.get_position(verbosity=0)

                if value > prev_value:
                    # Keep going in this direction...
                    pass
                else:
                    # Switch directions!
                    direction *= -1
                    step_size *= 0.5

        elif target == "min":

            if verbosity >= 5:
                print("Performing search on axis '{}' target is 'min'".format(self.name))

            direction = +1 * polarity

            while step_size >= min_step:
                if verbosity >= 4:
                    print("        move {} by {} × {}".format(self.name, direction, step_size))

                # pos = yield from bps.rd(motor)
                yield from bps.mvr(motor, direction * step_size)
                # self.move_relative(move_amount=direction*step_size, verbosity=verbosity-2)

                prev_value = value
                yield from bps.trigger_and_read([detector, *motors_for_table])
                value = detector.read()[value_name]["value"]
                if verbosity >= 3:
                    print(
                        "      {} = {:.3f} {}; value : {}".format(
                            self.name, self.get_position(verbosity=0), self.units, value
                        )
                    )

                if value < prev_value:
                    # Keep going in this direction...
                    pass
                else:
                    # Switch directions!
                    direction *= -1
                    step_size *= 0.5

        else:

            target_rel = target
            target = target_rel * intensity

            if verbosity >= 5:
                print(
                    "Performing search on axis '{}' target {} × {} = {}".format(
                        self.name, target_rel, intensity, target
                    )
                )
            if verbosity >= 4:
                print("      value : {} ({:.1f}%)".format(value, 100.0 * value / intensity))

            # Determine initial motion direction
            if value > target:
                direction = -1 * polarity
            else:
                direction = +1 * polarity

            while step_size >= min_step:

                if verbosity >= 4:
                    print("        move {} by {} × {}".format(self.name, direction, step_size))

                yield from bps.mvr(motor, direction * step_size)

                yield from bps.trigger_and_read([detector, *motors_for_table])
                value = detector.read()[value_name]["value"]
                # if verbosity>=3:
                #    print("      {} = {:.3f} {}; value : {} ({:.1f}%)".format(self.name, self.get_position(verbosity=0), self.units, value, 100.0*value/intensity))

                # Determine direction
                if value > target:
                    new_direction = -1.0 * polarity
                else:
                    new_direction = +1.0 * polarity

                if abs(direction - new_direction) < 1e-4:
                    # Same direction as we've been going...
                    # ...keep moving this way
                    pass
                else:
                    # Switch directions!
                    direction *= -1
                    step_size *= 0.5

        # bec.enable_table()

    # def calc_lookuptable(self,target_x):
    #     #make a look up table for 

    #     start_x = self.start_x
    #     start_y = self.start_y
    #     start_th = self.start_th
        
    #     end_x = self.end_x
    #     end_y = self.end_y
    #     end_th = self.end_th
       
    #     target_y = (target_x-end_x)/(start_x-end_x)*(start_y-end_y)+end_y
    #     target_th = (target_x-end_x)/(start_x-end_x)*(start_th-end_th)+end_th

    #     return target_x, target_y, target_th

    def run_initial_alignment(self,start_x=0, end_x=22, direct_beam_int=None):
        #make a look up table for 


        yield from bps.mv(smx, start_x)
        yield from self.align_crazy_v3_plan(direct_beam_int=direct_beam_int)

        start_x = smx.position
        start_y = smy.position
        start_th = sth.position

        yield from bps.mv(smx, end_x)
        yield from self.align_crazy_v3_plan(direct_beam_int=direct_beam_int)


        end_x = smx.position
        end_y = smy.position
        end_th = sth.position

        self.start_x = start_x
        self.start_y = start_y
        self.start_th = start_th
        self.end_x = end_x
        self.end_y = end_y
        self.end_th = end_th
        # start_x = self.start_x
        # start_y = self.start_y
        # start_th = self.start_th
        
        # end_x = self.end_x
        # end_y = self.end_y
        # end_th = self.end_th
        # return self.calc_lookuptable()

    def calc_lookuptable(self, target_x):

        start_x = self.start_x
        start_y = self.start_y
        start_th = self.start_th
        
        end_x = self.end_x
        end_y = self.end_y
        end_th = self.end_th

        target_y = (target_x-end_x)/(start_x-end_x)*(start_y-end_y)+end_y
        target_th = (target_x-end_x)/(start_x-end_x)*(start_th-end_th)+end_th

        return target_x, target_y, target_th


    def align_lookup(self, target_x, direct_beam_int=None):
        
        xpos, ypos, thpos = self.calc_lookuptable(target_x)
        
        #move to the position in lookup table
        yield from bps.mv(smx, xpos, smy, ypos, sth, thpos)
        
        yield from self.align_crazy_v3_plan(direct_beam_int=direct_beam_int)




    def search_plan2(
        self,
        motor=smy,
        step_size=0.2,
        min_step=0.05,
        intensity=None,
        target=0.5,
        detector=None,
        detector_suffix=None,
        polarity=-1,
        verbosity=3,
    ):
        """Moves this axis, searching for a target value.

        Parameters
        ----------
        step_size : float
            The initial step size when moving the axis
        min_step : float
            The final (minimum) step size to try
        intensity : float
            The expected full-beam intensity readout
        target : 0.0 to 1.0
            The target ratio of full-beam intensity; 0.5 searches for half-max.
            The target can also be 'max' to find a local maximum.
        detector, detector_suffix
            The beamline detector (and suffix, such as '_stats4_total') to trigger to measure intensity
        polarity : +1 or -1
            Positive motion assumes, e.g. a step-height 'up' (as the axis goes more positive)
        """

        if detector is None:
            # detector = gs.DETS[0]
            detector = get_beamline().detector[0]

        @bpp.stage_decorator([detector])
        @bpp.run_decorator(md={})
        @bpp.finalize_decorator(final_plan=shutter_off)
        def inner_search(group=None):
            nonlocal intensity, target, step_size

            if not get_beamline().beam.is_on():
                print("WARNING: Experimental shutter is not open.")

            if intensity is None:
                intensity = RE.md["beam_intensity_expected"]

            if group:
                yield from bps.wait(group)

            yield from shutter_on()

            yield from self.search_stub2(
                motor=motor,
                step_size=step_size,
                min_step=min_step,
                intensity=intensity,
                target=target,
                detector=detector,
                detector_suffix=detector_suffix,
                polarity=polarity,
                verbosity=verbosity,
            )

        group_name = "setup_aligment"
        yield from bps.abs_set(bsx, cms.bsx_pos + 3, group=group_name)
        # beam.setTransmission(1e-6)

        yield from inner_search(group=group_name)

        yield from bps.abs_set(bsx, cms.bsx_pos, group=group_name)
        yield from bps.wait(group_name)

    def measureAutonomous(
        self,
        runno=0,
        exposure_time=10,
        extra=None,
        max_measurements=600000,
        prefix="measureAutonomous > ",
        verbosity=3,
        **md,
    ):
        """Measure points in a loop, relying on an external queue to specify what
        position to actually measure. If the 'position' is not (x,y) sample coordinates,
        then you will have to add code to do the appropriate coordinate conversion,
        or trigger the right beamline motors/components."""

        for i in range(runno, max_measurements):

            if verbosity >= 3:
                print("{}Waiting for AE command on queue...".format(prefix))

            # forceload_repeat = 0
            # if forceLoad == True:
            #     commands = measure_queue.get()
            #     forceload_repeat = 1
            # elif forceLoad == False or forceload_repeat == 1:
            commands = measure_queue.get()  # Get measurement command from queue
            num_to_measure = sum([1.0 for command in commands if command["measured"] is False])

            if verbosity >= 3:
                # print('{}Received command to measure {} points'.format(num_to_measure))
                print("{}Received command to measure {} points".format(prefix, num_to_measure))

            imeasure = 0
            for icommand, command in enumerate(commands):
                if verbosity >= 5:
                    print("{}Considering point {}/{}".format(prefix, icommand, len(commands)))

                if not command["measured"]:
                    imeasure += 1
                    if verbosity >= 3:
                        print("{}Measuring point {}/{}".format(prefix, imeasure, num_to_measure))

                    start_time = time.time()

                    ########################################
                    # Move to point
                    ########################################
                    # Here you should define the beamline changes needed to go
                    # to the desired position. (You shouldn't in general need
                    # to change code outside of this block.)
                    ########################################

                    # convert x_pos, yy_pos of stage to x_position, y_position of sample
                    # yy_pos = -1*command['position']['x_position']
                    # x_pos = command['position']['y_position']

                    # [x_pos, yy_pos] = command['position']
                    Ti, Tm, Tf = command["position"]
                    # command['position_gpcam'] = command['position']

                    # align the sample
                    ss = input("Change the sample: (it has to be y or yes)")
                    while ss != "yes" and ss != "y":
                        ss = input("Change the sample: (it has to be y or yes)")
                    pta.laserOff()
                    smx.move(-20)

                    ss = input("Is the sample ready? (it has to be y or yes)")
                    while ss != "yes" and ss != "y":
                        ss = input("Is the sample ready? (it has to be y or yes)")

                    self.xo()
                    self.align()
                    cms.modeMeasurement()
                    self.xabs(0)
                    self.thabs(0.12)

                    # offset the sample for high T
                    # self.yr(xxx)
                    # self.thabs(xxx)

                    input("ready to go?")
                    # now use zmq to set PTA
                    # CustomQueue.push the message
                    BS.publish([Ti, Tm, Tf])
                    # check it's working
                    # if verbosity>=3:
                    #     print('{}Driving to point {}/{}; (x,yy) = ({:.3f}, {:.3f})'.format(prefix, imeasure, num_to_measure, x_pos, yy_pos))

                    self.name = (
                        self.name_o + "_Ti{:.1f}".format(Ti) + "_Tm{:.1f}".format(Tm) + "_Tf{:.1f}".format(Tf)
                    )

                    # continuous measurement for 5 points

                    wait_time_list = [5, 0, 0, 35, 85]
                    self.reset_clock()
                    # continuous measurement for 6 points
                    for ii in range(5):
                        self.xabs(ii * 0.2)
                        if ii >= 2:
                            self.crazy_th()
                            cms.modeMeasurement()
                        time.sleep(wait_time_list[ii])

                        self.measureIncidentAngle(0.12, exposure_time=exposure_time, **md)

                    # wait until the T back to RT, move to the last fresh position
                    while self.clock() < 310:
                        time.sleep(5)
                    self.xabs(0.2 * 6)
                    self.crazy_th()
                    cms.modeMeasurement()
                    self.measureIncidentAngle(0.12, exposure_time=exposure_time, extra="FINAL", **md)

                    header = db[-1]  # The most recent measurement
                    # command['filename'] = '{}'.format(header.start['filename'][:-1])
                    command["filename"] = "{}".format(header.start["filename"])

                    ########################################
                    # md['anneal_time'] = self.anneal_time
                    # md['preanneal_time'] = self.preanneal_time

                    cost_time = time.time() - start_time

                    command["cost"] = cost_time

                    command["measured"] = True
                    command["analyzed"] = False

            measure_queue.publish(commands)  # Send results for analysis

    def measureAutonomous_3samples(
        self,
        runNo=0,
        exposure_time=10,
        extra=None,
        max_measurements=600000,
        samplePosNo=0,
        prefix="measureAutonomous > ",
        verbosity=3,
        align=False,
        **md,
    ):
        """Measure points in a loop, relying on an external queue to specify what
        position to actually measure. If the 'position' is not (x,y) sample coordinates,
        then you will have to add code to do the appropriate coordinate conversion,
        or trigger the right beamline motors/components."""

        samplePosNo_list = np.arange(0, 3)
        samplePos = [123, 143.5, 163]
        # samplePos = self.smxPos
        # samplePosOffset = [pos1, pos2, pos3]
        samplePosNo = samplePosNo

        for i in range(runNo, max_measurements):

            if verbosity >= 3:
                print("{}Waiting for AE command on queue...".format(prefix))

            commands = measure_queue.get()  # Get measurement command from queue
            num_to_measure = sum([1.0 for command in commands if command["measured"] is False])

            if verbosity >= 3:
                # print('{}Received command to measure {} points'.format(num_to_measure))
                print("{}Received command to measure {} points".format(prefix, num_to_measure))

            imeasure = 0
            for icommand, command in enumerate(commands):
                if verbosity >= 5:
                    print("{}Considering point {}/{}".format(prefix, icommand, len(commands)))

                if not command["measured"]:
                    imeasure += 1
                    if verbosity >= 3:
                        print("{}Measuring point {}/{}".format(prefix, imeasure, num_to_measure))

                    start_time = time.time()

                    # [x_pos, yy_pos] = command['position']
                    Ti, Tm, Tf = command["position"]
                    print("The next Temperatures are: Ti {:.1f}, Tm {:.1f}, Tf {:.1f}".format(Ti, Tm, Tf))
                    # command['position_gpcam'] = command['position']
                    self.name = (
                        self.name_o
                        + "_Ti{:.1f}".format(Ti)
                        + "_Tm{:.1f}".format(Tm)
                        + "_Tf{:.1f}".format(Tf)
                        + "_run{}".format(i)
                    )

                    # align the sample
                    if samplePosNo == 0:

                        ss = input("Change the sample: (it has to be y or yes)")
                        while ss != "yes" and ss != "y":
                            ss = input("Change the sample: (it has to be y or yes)")
                        pta.laserOff()
                        smx.move(-20)

                        ss = input("Are the samples ready? (it has to be y or yes)")
                        while ss != "yes" and ss != "y":
                            ss = input("Are the samples ready? (it has to be y or yes)")
                    print("Running Postion {}".format(samplePosNo))
                    smx.move(samplePos[samplePosNo])
                    sth.move(0)
                    self.setOrigin(["th"])

                    while smx.moving == True:
                        time.sleep(1)
                    beam.setSize(0.2, 0.05)
                    self.align()
                    beam.setSize(0.2, 0.2)
                    self.setOrigin(["x"])
                    input("ready to go?")
                    # now use zmq to set PTA
                    # CustomQueue.push the message
                    BS.publish([Ti, Tm, Tf])
                    # check it's working

                    # continuous measurement for 5 points

                    wait_time_list = [5, 0, 0, 50, 100]

                    self.reset_clock()
                    for ii in range(5):
                        self.xabs(ii * 0.2)
                        if ii >= 2 and align:
                            self.crazy_th()
                        time.sleep(wait_time_list[ii])
                        cms.modeMeasurement()
                        self.measureIncidentAngle(0.12, exposure_time=exposure_time)

                        # wait until the T back to RT, move to the last fresh position
                    while self.clock() < 310:
                        time.sleep(5)
                    self.xabs(0.2 * 5)
                    if align:
                        self.crazy_th()
                    else:
                        beam.setSize(0.2, 0.05)
                        self.align()
                        beam.setSize(0.2, 0.2)
                    cms.modeMeasurement()
                    self.measureIncidentAngles(
                        [0.08, 0.1, 0.12, 0.15, 0.18, 0.2, 0.25], exposure_time=exposure_time, extra="FINAL"
                    )

                    header = db[-3]  # The most recent measurement
                    # command['filename'] = '{}'.format(header.start['filename'][:-1])
                    command["filename"] = "{}".format(header.start["filename"])

                    ########################################
                    # md['anneal_time'] = self.anneal_time
                    # md['preanneal_time'] = self.preanneal_time

                    cost_time = time.time() - start_time

                    command["cost"] = cost_time

                    command["measured"] = True
                    command["analyzed"] = False

                    samplePosNo += 1
                    if samplePosNo > 2:
                        samplePosNo = 0

            measure_queue.publish(commands)  # Send results for analysis

    def measureManual(self, Ti, Tm, Tf, step=0, exposure_time=10, align=False):
        if step < 1:
            ss = input("Sample pos number? ")
            if ss == "0":
                smx.move(self.smxPos[0])
            elif ss == "1":
                smx.move(self.smxPos[1])
            elif ss == "2":
                smx.move(self.smxPos[2])
            else:
                print("Wrong Number. It has to be 0 or 1 or 2.")
                return

            while smx.moving == True:
                time.sleep(1)

        if step < 5:
            self.setOrigin(["x"])
            self.gotoOrigin()
            beam.setSize(0.2, 0.05)
            self.align()
            beam.setSize(0.2, 0.2)

        if step < 10:
            ss = input("Ready for annealing? (it has to be y or yes)")
            while ss != "yes" and ss != "y":
                ss = input("Ready for annealing? (it has to be y or yes)")

            BS.publish([Ti, Tm, Tf])
            if align:
                wait_time_list = [5, 0, 0, 35, 85]
            else:
                wait_time_list = [5, 0, 0, 50, 100]
            self.reset_clock()
            # continuous measurement for 6 points
            for ii in range(5):
                self.xabs(ii * 0.2)
                if ii >= 2 and align:
                    self.crazy_th()
                time.sleep(wait_time_list[ii])
                cms.modeMeasurement()
                self.measureIncidentAngle(0.12, exposure_time=exposure_time)

                # wait until the T back to RT, move to the last fresh position
            while self.clock() < 310:
                time.sleep(5)
            self.xabs(0.2 * 5)
            if align:
                self.crazy_th()
            else:
                beam.setSize(0.2, 0.05)
                self.align()
                beam.setSize(0.2, 0.2)
            cms.modeMeasurement()
            self.measureIncidentAngles(
                [0.08, 0.1, 0.12, 0.15, 0.18, 0.2, 0.25], exposure_time=exposure_time, extra="FINAL"
            )

    def alignment_set(self):
        self.start_x = 0
        self.end_x = 22
        self.start_y = 39.2196
        self.end_y = 40.037
        self.start_th = 1.2
        self.end_th = 1.198
        cms.direct_beam_int = 22245
        yield from bps.null()



# cms.SAXS.setCalibration([758, 1680-607], 5.0, [-65, -73]) # 13.5 keV
# cms.SAXS.setCalibration([754, 1076], 5.03, [-65, -73])  #20190320, 13.5 keV
# cms.SAXS.setCalibration([754, 1075], 5.03, [-65, -73])  #20201021, 13.5 keV
# cms.SAXS.setCalibration([754, 1077], 5.03, [-65, -73])  #20210208, 13.5 keV
# cms.SAXS.setCalibration([761, 1075], 5.03, [-65, -73])  #20210716, 13.5 keV
cms.SAXS.setCalibration([756, 1079], 5.83, [-65, -73])  # 20201021, 13.5 keV


# RE.md['experiment_group'] = 'MNoack'
RE.md["experiment_group"] = "KYager"
# RE.md['experiment_alias_directory'] = '/nsls2/xf11bm/data/2020_3/MNoack/Exp1/'
RE.md["experiment_alias_directory"] = "/nsls2/data/cms/legacy/xf11bm/data/2023_1/beamline/PTA2/"
RE.md["experiment_user"] = "TBD"
RE.md["experiment_type"] = "SAXS"
RE.md["experiment_project"] = "TBD"


def fake_coordinated_motion(mtr1, target1, mtr2, target2, N=1000):
    start1 = yield from bps.rd(mtr1)
    start2 = yield from bps.rd(mtr2)
    step1 = (target1 - start1) / N
    step2 = (target2 - start2) / N
    for j in range(int(N)):
        yield from bps.mv(mtr1, start1 + j * step1, mtr2, start2 + j * step2)


def fake_coordinated_motionr(mtr1, mtr2, delta, step=0.1):
    
    real_step= step*abs(delta)/delta

    for j in range(int(abs(delta)/ step)):
        
        yield from bps.mvr(mtr1, real_step, mtr2, real_step)


def changeSamplesa():
    smx.move(-20)


# the strips in the same materials system should have the same file name
# the location with l_pos=0, o_pos=0 for the first strip should be mared as 0
# the location with l_pos=0, o_pos=i for the second sample needs to be aligned in the sequence of the first strip.


# # Connect to ZeroMQ (zmq)
# try:
#     BS
# except NameError:
#     ##queue_PATH='../'
#     queue_PATH='/nsls2/data/cms/legacy/xf11bm/data/2022_2/SRussell/'
#     queue_PATH in sys.path or sys.path.append(queue_PATH)
#     from CustomQueuePTA import BSQueue
#     BS = BSQueue()

# ## Connect to S3
# try:
#     measure_queue
# except NameError:
#     ##queue_PATH='../'
#     queue_PATH='/nsls2/data/cms/legacy/xf11bm/data/2022_2/SRussell/'
#     queue_PATH in sys.path or sys.path.append(queue_PATH)
#     from CustomS3 import Queue_measure
#     measure_queue = Queue_measure()

"""
The edge of the diving board in x

In [237]: wsam()
smx = 81.2
smy = 16.13191875
sth = 1.120000000000001

the other edge of x (clamping spot)

smx = 106.1

============
align bare Si wafer at smx=106

In [298]: wsam()
smx = 106.0
smy = 16.19016875
sth = 0.995000000000001

align bare Si wafer at smx=81.1

In [298]: wsam()

In [365]: wsam()
smx = 93.55
smy = 15.7754375
sth = 0.9976562500000004

1.89deg offset in schi 


============
y scan at the aligned position to verify the stats2/stat4.
+-----------+------------+------------+-------------------+------------------------+------------------------+------------------------+------------------------+
|   seq_num |       time |        smy | smy_user_setpoint | pilatus2M_stats1_total | pilatus2M_stats2_total | pilatus2M_stats3_total | pilatus2M_stats4_total |
+-----------+------------+------------+-------------------+------------------------+------------------------+------------------------+------------------------+
|         1 | 12:02:13.5 |    15.5819 |           15.5819 |                      1 |                    -89 |                      0 |                  20095 |
|         2 | 12:02:15.4 |    15.6019 |           15.6019 |                      1 |                    -90 |                      0 |                  20221 |                 
|         3 | 12:02:17.3 |    15.6219 |           15.6219 |                      0 |                    -90 |                      0 |                  20026 |                 
|         4 | 12:02:19.1 |    15.6419 |           15.6419 |                      0 |                    -90 |                      0 |                  20337 |                 
|         5 | 12:02:21.1 |    15.6619 |           15.6619 |                      1 |                    -90 |                      0 |                  20427 |                 
|         6 | 12:02:22.9 |    15.6819 |           15.6819 |                      2 |                    -87 |                      0 |                  20053 |                 
|         7 | 12:02:24.7 |    15.7019 |           15.7019 |                      2 |                    -84 |                      3 |                  20383 |                 
|         8 | 12:02:26.5 |    15.7219 |           15.7219 |                     19 |                    -51 |                      3 |                  20052 |                 
|         9 | 12:02:28.4 |    15.7419 |           15.7419 |                    186 |                    328 |                     75 |                  19790 |                 
|        10 | 12:02:30.4 |    15.7619 |           15.7619 |                   1450 |                   3553 |                    718 |                  16054 |                 
|        11 | 12:02:32.3 |    15.7819 |           15.7819 |                   3974 |                   8884 |                   2697 |                   7412 |                 
|        12 | 12:02:34.2 |    15.8019 |           15.8019 |                   2783 |                   5922 |                   1889 |                   1563 |                 
|        13 | 12:02:36.1 |    15.8219 |           15.8219 |                    705 |                   1368 |                    558 |                    114 |                 
|        14 | 12:02:38.0 |    15.8419 |           15.8419 |                     65 |                     39 |                     61 |                      2 |                 
|        15 | 12:02:39.9 |    15.8619 |           15.8619 |                     15 |                    -66 |                     10 |                      0 |                 
|        16 | 12:02:41.8 |    15.8819 |           15.8819 |                      7 |                    -70 |                      0 |                      0 |                 
|        17 | 12:02:43.5 |    15.9019 |           15.9019 |                      0 |                    -89 |                      0 |                      0 |                 
|        18 | 12:02:45.5 |    15.9219 |           15.9219 |                      0 |                    -90 |                      0 |                      0 |                 
|        19 | 12:02:47.4 |    15.9419 |           15.9419 |                      0 |                    -90 |                      0 |                      0 |                 
|        20 | 12:02:49.4 |    15.9619 |           15.9619 |                      0 |                    -90 |                      0 |                      0 |                 
|        21 | 12:02:51.3 |    15.9819 |           15.9819 |                      0 |                    -90 |                      0 |                      0 |                 



+-----------+------------+------------+-------------------+------------------------+------------------------+------------------------+------------------------+


laser power <24

0, 6, 12, 18, 24

heater @50C @100C.


#laser @ the edge of Si wafer 
In [850]: wsam()
smx = 86.2495
smy = 15.778553125
sth = 1.0212500000000002

#2mm offset from laser @ the edge
In [854]: wsam()
smx = 88.2495
smy = 15.778553125
sth = 1.0212500000000002
In [856]: laserx.position
Out[856]: 0.0

FIX the offfset between smx and laserx as 88.25


#alignment scan @50C


# smy.mov(15.5)
RE(fake_coordinated_motion(smx, 82, laserx, -6.25, N=120))
#set power
for power in np.arange(0, 24.1, 6):
    
    pta.setLaserPower(power)
    #set x position
    if power==0:
        pta.laserOff()
    else:
        pta.laserOn()
        time.sleep(30)
    # for xpos in np.arrange(0, 24.1, 6):
    for xpos in range(5):
        smy.move(15.5)
        sam.align()
        if xpos<4:
            RE(fake_coordinated_motionr(smx, laserx, delta=6))
    
    pta.laserOff()

    RE(fake_coordinated_motion(smx, 82, laserx, -6.25, N=120))


pta.setLaserPower(power)


"""


# camonitor -S XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire

# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:18:45.322852 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:18:45.406499 Done  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:18:45.332771 23333  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:18:45.332854 10  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:18:45.562211 Done  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:10.780798 0  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:11.563899 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:11.563945 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:11.563973 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:11.632702 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:11.937539 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:11.953721 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:11.953786 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:11.963590 23334  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:11.963779 1  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:12.055200 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:12.278970 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:12.279087 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:12.279116 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:12.347917 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:12.653034 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:12.666444 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:12.666491 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:12.674896 23335  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:12.674989 2  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:12.806157 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:12.929527 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:12.929575 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:12.929593 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:12.998175 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:13.301574 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:13.316497 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:13.316522 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:13.325644 23336  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:13.325764 3  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:13.452583 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:13.567623 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:13.567675 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:13.567699 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:13.636079 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:14.153505 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:14.168763 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:14.168805 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:14.178550 23337  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:14.178706 4  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:14.252755 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:14.380508 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:14.380533 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:14.380543 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:14.449284 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:14.756593 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:14.772070 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:14.772109 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:14.781253 23338  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:14.781299 5  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:14.859708 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:14.982903 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:14.982950 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:14.982985 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:15.051642 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:15.661136 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:15.676297 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:15.676319 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:15.685042 23339  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:15.685210 6  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:15.816051 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:15.932728 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:15.932767 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:15.932785 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:16.001699 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:16.305201 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:16.319495 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:16.319542 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:16.329081 23340  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:16.329159 7  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:16.408806 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:16.533301 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:16.533346 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:16.533364 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:16.602075 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:17.143514 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:17.157371 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:17.157411 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:17.167105 23341  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:17.167203 8  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:17.245183 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:17.369047 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:17.369103 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:17.369129 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:17.437877 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:17.743062 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:17.756966 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:17.757005 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:17.766715 23342  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:17.766777 9  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:17.919766 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:18.042324 Acquiring data  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:18.042377 Acquire  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:18.042401 Acquiring STATE MINOR
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:18.111180 Waiting for 7OK response  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:18.673632 Reading image file /ramdisk/current_0000.tiff  
# XF:11BMB-ES{Det:PIL2M}:cam1:Acquire 2023-03-20 15:25:18.688266 Done  
# XF:11BMB-ES{Det:PIL2M}:cam1:StatusMessage_RBV 2023-03-20 15:25:18.688318 Waiting for acquire command  
# XF:11BMB-ES{Det:PIL2M}:Trans1:ArrayCounter_RBV 2023-03-20 15:25:18.697376 23343  
# XF:11BMB-ES{Det:PIL2M}:TIFF1:ArrayCounter_RBV 2023-03-20 15:25:18.697448 10  
# XF:11BMB-ES{Det:PIL2M}:cam1:AcquireBusy 2023-03-20 15:25:18.825073 Done  


def test_plan(detector=None):

    if detector is None:
        detector=pilatus2M
        # detector = get_beamline().detector[0]

    motors_for_table = [smx, smy, sth]

    @bpp.stage_decorator([detector])
    @bpp.run_decorator(md={})
    @bpp.finalize_decorator(final_plan=shutter_off)
    def inner_plan(group=None):

        if group:
            yield from bps.wait(group)

        for n in range(10):
            t0 = time.time()
            # yield from bps.trigger_and_read([detector, *motors_for_table])
            yield from bps.trigger_and_read([detector])
            print(f"Detection time: {time.time() - t0}")
            # yield from bps.sleep(.1)

    group_name = "setup_aligment"
    yield from bps.abs_set(bsx, cms.bsx_pos + 3, group=group_name)
    beam.setTransmission(1e-6)

    yield from inner_plan(group=group_name)

    yield from bps.abs_set(bsx, cms.bsx_pos, group=group_name)
    yield from bps.wait(group_name)

sam = Sample('test')

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



'''
Notes at Mar 22, 2023, the 4th day at CMS

#quick alignment for Si wafer

#check the beam position in USB cam1
#

#the start edge 
In [165]: wsam()
smx = 0.0
smy = 39.23188125
sth = 1.086875000000001

In [166]: sam.pos()
test.th = 0.157 deg (origin = 0.930)
test.x = 0.000 mm (origin = 0.000)
test.y = -0.028 mm (origin = 39.260)
Out[166]: {'th': 0.15687500000000087, 'x': 0.0, 'y': -0.028118749999997306}


#the end edge

In [173]: wsam()
smx = 22.0005
smy = 40.072500000000005
sth = 1.0865624999999994

In [174]: sam.pos()
test.th = 0.157 deg (origin = 0.930)
test.x = 22.000 mm (origin = 0.000)
test.y = 0.813 mm (origin = 39.260)

'''


'''
Notes March 22, 2023
17:18 KY loaded "real scientific sample
name='GD4-113-4-2nd_Tb50C'

laserPower = 1.4V (4.6w)

smx = 3.7 ; laserx = 11.5 # IR laser is hitting edge of sample
smx = 8.7 ; laserx = 11.5 # IR laser is hitting edge of sample


Notes March 23, 2023
10:40 KY loaded "real scientific sample 2
name='GD4-113-1_Tb50C'
laserPower = 1.1V (2.3w)

there is unexpected ~10s delay on every point during alignment. the bug is gone after restart. 

10:40 KY loaded "real scientific sample 3
name='GD4-113-1-2nd_Tb50C'
laserPower = 1.1V (2.3w)

In [7]: sam.start_x
Out[7]: 0.0

In [8]: sam.start_y
Out[8]: 39.218578125

In [9]: sam.start_th
Out[9]: 1.2006250000000005

In [10]: sam.end_x
Out[10]: 22.0005

In [11]: sam.end_y
Out[11]: 40.037328125

In [12]: sam.end_th
Out[12]: 1.1979687500000011


#####protocol
=============bsui==================
#change sample --mov smx -100 and mount the fresh sample
RE(sam.run_initial_alignment(start_x=0, end_x=22))  #align samples at smx=0 and 22 and calculate the lookup table for smy and th
#print out sam.start_x/y/th and smx.end_x/y/th and HARD code them in sam.alignment_set()
=============bsui part is done. exit bsui=============================
#restart the env in Queue monitor
#pre-load 'agent_alignemnt_set' and 'agent_laser_on'
==>>ws2, agent.start(True)
==>>ws1, start the queue
'''

def fake_fly(mtr, start, stop, exp_time):
    det = pilatus2M

    det.tiff.kind = 'omitted'
    det.tiff.disable_on_stage()
    det.stats4.total.kind='hinted'
   

    @bpp.stage_decorator([det])
    @bpp.monitor_during_decorator([det.stats4.total, mtr])
    @bpp.run_decorator()
    @bpp.finalize_decorator(final_plan=shutter_off)
    def inner():
        yield from shutter_on()

        yield from bps.trigger(det, group='fake_fly')
        yield from bps.abs_set(mtr, stop, group='fake_fly')
        yield from bps.wait(group='fake_fly')

    @bpp.reset_positions_decorator([det.cam.num_images, det.cam.acquire_time, det.cam.acquire_period])
    def inner2():
        yield from bps.mv(det.cam.acquire_time, exp_time)
        yield from bps.mv(det.cam.acquire_period, exp_time +.05)
        total_time = np.abs(stop - start)

        yield from bps.mv(det.cam.num_images, num)
        yield from bps.mv(mtr, start)
        yield from inner()

    group_name = "setup_aligment"
    yield from bps.abs_set(bsx, cms.bsx_pos + 3, group=group_name)
    yield from bps.wait(group=group_name)

    beam.setTransmission(1e-6)

    yield from inner2()

    yield from bps.abs_set(bsx, cms.bsx_pos, group=group_name)
    yield from bps.wait(group_name)


def fake_fly2_test(mtr, start, stop, step, exp_time):

    # motors: smy (+/- 2), sth (+/- 1)
    det = pilatus2M

    # It takes 0.4 to 0.7 s longer to complete motion, so let's add 1 s for now
    #   It should be computed/estimated more accurately
    num = int(np.abs(stop - start) / step)
    total_time = num * exp_time
    velocity = np.abs(stop - start) / total_time

    det.tiff.kind = 'omitted'
    det.tiff.disable_on_stage()
    det.stats4.total.kind='hinted'

    frame_numbers = []
    frame_timestamps = []
    frame_mtr_pos = []
    frame_roi2_int = []
    frame_roi3_int = []
    frame_roi4_int = []

    def accumulate(value, old_value, timestamp, **kwargs):
        frame_numbers.append(value)
        frame_timestamps.append(timestamp)
        frame_mtr_pos.append(mtr.position)

    def accumulate_roi2(value, old_value, timestamp, **kwargs):
        roi2_int = pilatus2M.stats2.total.get()
        frame_roi2_int.append(roi2_int)

    def accumulate_roi3(value, old_value, timestamp, **kwargs):
        roi3_int = pilatus2M.stats3.total.get()
        frame_roi3_int.append(roi3_int)

    def accumulate_roi4(value, old_value, timestamp, **kwargs):
        roi4_int = pilatus2M.stats4.total.get()
        frame_roi4_int.append(roi4_int)

    @bpp.stage_decorator([det])
    def inner():
        cid = pilatus2M.cam.array_counter.subscribe(accumulate)
        cid2 = pilatus2M.stats2.array_counter.subscribe(accumulate_roi2)
        cid3 = pilatus2M.stats3.array_counter.subscribe(accumulate_roi3)
        cid4 = pilatus2M.stats4.array_counter.subscribe(accumulate_roi4)
        try:
            yield from bps.trigger(det, group='fake_fly')
            yield from bps.abs_set(mtr, stop, group='fake_fly')
            yield from bps.wait(group='fake_fly')
        finally:
            pilatus2M.cam.array_counter.unsubscribe(cid)
            pilatus2M.stats2.array_counter.unsubscribe(cid2)
            pilatus2M.stats3.array_counter.unsubscribe(cid3)
            pilatus2M.stats4.array_counter.unsubscribe(cid4)

    @bpp.reset_positions_decorator([det.cam.num_images, det.cam.acquire_time, det.cam.acquire_period,
                                    mtr.velocity])
    def inner2():
        yield from bps.abs_set(mtr, start, wait=True)
        yield from bps.mv(mtr.velocity, velocity)

        print(f"Number of acquired images: {num}. Exposure time: {exp_time}")
        
        yield from bps.mv(det.cam.acquire_time, exp_time - 0.005)
        yield from bps.mv(det.cam.acquire_period, exp_time)

        yield from bps.mv(det.cam.num_images, num)
        yield from inner()

    yield from inner2()

    def trim_list(v, num):
        n_first = max(len(v) - num, 0)
        return v[n_first:]

    frame_numbers = trim_list(frame_numbers, num)
    frame_timestamps = trim_list(frame_timestamps, num)
    frame_mtr_pos = trim_list(frame_mtr_pos, num)
    frame_roi2_int = trim_list(frame_roi2_int, num)
    frame_roi3_int = trim_list(frame_roi3_int, num)
    frame_roi4_int = trim_list(frame_roi4_int, num)
    
    print(f"frame_numbers = {frame_numbers}")
    print(f"frame_timestamps = {frame_timestamps}")
    print(f"mtr_pos = {frame_mtr_pos}")
    print(f"roi2 = {frame_roi2_int}")
    print(f"roi3 = {frame_roi3_int}")
    print(f"roi4 = {frame_roi4_int}")

    return frame_mtr_pos, frame_roi2_int, frame_roi3_int, frame_roi4_int


def agent_feedback_plan(sample_x, md=None):
    md = md or {}

    yield from sam.align_lookup(sample_x)
    print("ALIGN DONE")
    yield from cms.modeMeasurement_plan()
    print("IN MEASUE MODE")
    yield from sam.measure(1, **md)
    print("DONE")

def agent_bootstrap_alignment():
    yield from sam.run_initial_alignment()