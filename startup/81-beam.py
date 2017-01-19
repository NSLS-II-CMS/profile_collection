#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4




################################################################################
#  Code for querying and controlling beamline components that 'affect' the
# beam. (Beam energy, beam flux, etc.)
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################


# Notes
################################################################################
# verbosity=0 : Output nothing
# verbosity=1 : Output only final (minimal) result
# verbosity=2 : Output 'regular' amounts of information/data
# verbosity=3 : Output all useful information
# verbosity=4 : Output marginally useful things (e.g. essentially redundant/obvious things)
# verbosity=5 : Output everything (e.g. for testing)



# These imports are not necessary if part of the startup sequence.
# If this file is called separately, some of these may be needed.
#import numpy as np
#from epics import caget, caput
#from time import sleep

#from ophyd import EpicsMotor, Device, Component as Cpt
#from ophyd.commands import * # For mov, movr


class BeamlineElement(object):
    '''Defines a component of the beamline that (may) intersect the x-ray beam.'''
    
    def __init__(self, name, zposition, description="", pv=None, **args):
        
        self.name = name
        self.zposition = zposition
        self.description = description
        
        self.conversion_factor = 1
        
        self._pv_main = pv
        
        self.has_flux = True
        
        
    def state(self):
        """
        Returns the current state of the beamline element. Common states:
        out - Element is out of the way of the beam (and should not be blocking).
        in - Element is in the beam (but should not be blocking).
        block - Element is in the beam, and should be blocking the beam.
        undefined - Element is in an unexpected state.
        """
        
        return "out"

    
    def transmission(self, verbosity=0):
        """
        Returns the predicted transmission of this beamline element, based on 
        its current state.
        """
        
        tr_tot = 1.0
        
        if verbosity>=2:
            print('{:s} transmission = {:.6g}'.format(self.name, tr_tot))
        
        
        # Assume a generic beamline element doesn't block/perturb the beam
        return tr_tot
        
        
    def flux(self, verbosity=3):
        
        reading = self.reading(verbosity=0)
        flux = self.conversion_factor*reading # ph/s
        
        if verbosity>=2:
            print('flux = {:.4g} ph/s'.format(flux))
        
        return flux
    
    
        
        
class Shutter(BeamlineElement):
    
    # Example
    #          XF:11BMA-PPS{PSh}Enbl-Sts
    #  Status: XF:11BMA-PPS{PSh}Pos-Sts       0 for open, 1 for close
    #  Open:   XF:11BMA-PPS{PSh}Cmd:Opn-Cmd
    #  Close:  XF:11BMA-PPS{PSh}Cmd:Cls-Cmd
    
    def __init__(self, name, zposition, description="", pv=None, **args):
        
        super().__init__(name=name, zposition=zposition, description=description, pv=pv, **args)
        self.has_flux = False
        
    
    def state(self):
        """
        Returns the current state of the beamline element. Common states:
        out - Element is out of the way of the beam (and should not be blocking).
        in - Element is in the beam (but should not be blocking).
        block - Element is in the beam, and should be blocking the beam.
        undefined - Element is in an unexpected state.
        """
        
        state_n = caget(self._pv_main+'Pos-Sts')
        
        if state_n is 0:
            return "out"
        elif state_n is 1:
            return "block"
        else:
            return "undefined" 
        
        
    def open(self, verbosity=3):
        
        if verbosity>=3:
            print('Opening {:s}...'.format(self.name))
        
        # E.g. #XF:11BMB-VA{Slt:4-GV:1}Cmd:Opn-Cmd
        pv = self._pv_main + 'Cmd:Opn-Cmd'
        #caput(pv, 1) # TODO: Test this.
    
    def close(self, verbosity=3):
        
        if verbosity>=3:
            print('Closing {:s}...'.format(self.name))
            
        pv = self._pv_main + 'Cmd:Cls-Cmd'
        #caput(pv, 1) # TODO: Test this.

        



class GateValve(Shutter):
    
    # Example
    #  Status: XF:11BMB-VA{Slt:4-GV:1}Pos-Sts        1 for open, 0 for close
    #  Open:   XF:11BMB-VA{Slt:4-GV:1}Cmd:Opn-Cmd
    #  Close:  XF:11BMB-VA{Slt:4-GV:1}Cmd:Cls-Cmd
    
    
    def state(self):
        """
        Returns the current state of the beamline element. Common states:
        out - Element is out of the way of the beam (and should not be blocking).
        in - Element is in the beam (but should not be blocking).
        block - Element is in the beam, and should be blocking the beam.
        undefined - Element is in an unexpected state.
        """
        
        state_n = caget(self._pv_main+'Pos-Sts')
        
        if state_n is 1:
            return "out"
        elif state_n is 0:
            return "block"
        else:
            return "undefined"     
    


class ThreePoleWiggler(BeamlineElement):
    
    def __init__(self, name='3PW', zposition=0.0, description='Three-pole wiggler source of x-rays', **args):
        
        
        super().__init__(name=name, zposition=zposition, description=description, **args)
        
        # TODO: Find out the right conversion factor
        self.conversion_factor = 3e18/500.0 #(ph/s)/mA
        

    def state(self):
        """
        Returns the current state of the beamline element. Common states:
        out - Element is out of the way of the beam (and should not be blocking).
        in - Element is in the beam (but should not be blocking).
        block - Element is in the beam, and should be blocking the beam.
        undefined - Element is in an unexpected state.
        """
        
        position = caget('SR:C11-ID:G5{3PW:1}Mtr.RBV')
        
        # TODO: Instead use the 'inserted' flag?
        # caget('SR:C11-ID:G5{3PW:1}InsertedFlag')
        
        if abs(position-0)<3:
            return "in"
        
        elif abs(position - -189.0)<10:
            return "out"
        
        else:
            return "undefined"
        
        
    def reading(self, verbosity=3):
        
        if self.state() is 'in':
            
            ring_current = caget('SR:OPS-BI{DCCT:1}I:Real-I')
            if verbosity>=2:
                print('{:s} is inserted; ring current = {:.1f} mA'.format(self.name, ring_current))
                
            return ring_current
        
        else:
            if verbosity>=2:
                print('{:s} is not inserted.'.format(self.name))
                
            return 0
        

class Monitor(BeamlineElement):
    
    def quickReading(self, verbosity=3, delay=1.0):
        """
        Puts the diagnostic into the beam, takes a reading, and removes the
        diagnostic.
        """
        
        self.insert()
        sleep(delay)
        value = self.reading(verbosity=verbosity)
        
        self.retract()
        sleep(delay)
        
        return value
    
    
    
class DiagnosticScreen(Monitor):
    
    #XF:11BMB-BI{FS:4}Pos-Sts
    
    def __init__(self, name, zposition, description="", pv=None, epics_signal=None, **args):
        
        super().__init__(name=name, zposition=zposition, description=description, pv=pv, **args)
        self.epics_signal = epics_signal
        self.has_flux = False
        
    
    def state(self):
        """
        Returns the current state of the beamline element. Common states:
        out - Element is out of the way of the beam (and should not be blocking).
        in - Element is in the beam (but should not be blocking).
        block - Element is in the beam, and should be blocking the beam.
        undefined - Element is in an unexpected state.
        """
        
        state_n = caget(self._pv_main+'Pos-Sts')
        
        if state_n is 0:
            return "out"
        elif state_n is 1:
            return "block"
        else:
            return "undefined" 
            
    
    def insert(self, verbosity=3):
        
        if verbosity>=3:
            print('Inserting {:s}...'.format(self.name))
        
        # E.g. #XF:11BMB-VA{Slt:4-GV:1}Cmd:Opn-Cmd
        pv = self._pv_main + 'Cmd:In-Cmd'
        caput(pv, 1)
    
    def retract(self, verbosity=3):
        
        if verbosity>=3:
            print('Retracting {:s}...'.format(self.name))
            
        pv = self._pv_main + 'Cmd:Out-Cmd'
        caput(pv, 1)
        
        
    def reading(self, verbosity=3):
        
        value = self.epics_signal.stats1.total.value
        
        if self.state() is 'block':
            
            ring_current = caget('SR:OPS-BI{DCCT:1}I:Real-I')
            if verbosity>=2:
                print('{:s} is inserted; reading = {:.4g}'.format(self.name, value))
                
            return value
        
        else:
            if verbosity>=2:
                print('{:s} is not inserted.'.format(self.name))
                
            return 0
        
        
        
        
class PointDiode_CMS(Monitor):
    
    def __init__(self, name='bim6 point diode', zposition=58.3, description="Bar holding a point-diode, downstream of sample.", pv='XF:11BMB-BI{IM:2}EM180:Current1:MeanValue_RBV', epics_signal=None, **args):
        
        super().__init__(name=name, zposition=zposition, description=description, pv=pv, **args)
        self.has_flux = True
        
        if epics_signal is None:
            
            #bim6 = EpicsSignalROWait("XF:11BMB-BI{IM:2}EM180:Current1:MeanValue_RBV", wait_time=1, name='bim6')
            #bim6_integrating = EpicsSignalROIntegrate("XF:11BMB-BI{IM:2}EM180:Current1:MeanValue_RBV", wait_time=0.5, integrate_num=8, integrate_delay=0.1, name='bim6')
            
            self.epics_signal = bim6_integrating
            
        else:
            self.epics_signal = epics_signal
        
        
        # The beam (at the ion chamber) is roughly 0.50x0.50 mm.
        # If we slit down to 0.20x0.05 mm, we are capturing 0.4*0.25 = 0.1 of the beam.
        # bim6 reads 70000 cts (of course this depends on settings) when ion chamber reads 1.3e11 ph/s.
        # (settings: trans = 5e-4)
        # So conversion_factor is roughly:
        self.conversion_factor = 1.3e11*0.1/70000. # (ph/s)/cts
        
        self.in_position_x = 0.0
        self.in_position_y = 0.0

        self.out_position_x = 0.0
        self.out_position_y = -16.0
        
        self.position_tolerance = 0.1
        
    
    def state(self):
        """
        Returns the current state of the beamline element. Common states:
        out - Element is out of the way of the beam (and should not be blocking).
        in - Element is in the beam (but should not be blocking).
        block - Element is in the beam, and should be blocking the beam.
        undefined - Element is in an unexpected state.
        """
        
        position_x = DETx.user_readback.value
        position_y = DETy.user_readback.value
        
        if abs(position_x-self.out_position_x)<self.position_tolerance and abs(position_y-self.out_position_y)<self.position_tolerance: 
            return "out"
        if abs(position_x-self.in_position_x)<self.position_tolerance and abs(position_y-self.in_position_y)<self.position_tolerance: 
            return "block"
        else:
            return "undefined" 
                
    def insert(self, verbosity=3):
        
        if verbosity>=3:
            print('Inserting {:s}...'.format(self.name))
        
        mov( [DETx, DETy], [self.in_position_x, self.in_position_y] )
        
    
    def retract(self, verbosity=3):
        
        if verbosity>=3:
            print('Retracting {:s}...'.format(self.name))
            
        mov( [DETx, DETy], [self.out_position_x, self.out_position_y] )
        
        
    def reading(self, verbosity=3):
        
        value = self.epics_signal.read()[self.epics_signal.name]['value']
        
        if self.state() is 'block':
            
            if verbosity>=2:
                print('{:s} is inserted; reading = {:.4g}'.format(self.name, value))
                
            return value
        
        else:
            if verbosity>=2:
                print('{:s} is not inserted.'.format(self.name))
                
            return value
                
        

class IonChamber_CMS(Monitor):
    
    def __init__(self, name='bim3 ionchamber', zposition=49, description="Ion chamber (FMB Oxford I404) at start of endstation hutch", pv=None, beam=None, **args):
        
        super().__init__(name=name, zposition=zposition, description=description, pv=pv, **args)
        self.has_flux = True
        
        self.beam = beam
        
        # PVs
        import epics
        self.v1 = epics.PV('XF:11BMB-BI{IM:3}:IC1_MON')
        self.v2 = epics.PV('XF:11BMB-BI{IM:3}:IC2_MON')
        self.h1 = epics.PV('XF:11BMB-BI{IM:3}:IC3_MON')
        self.h2 = epics.PV('XF:11BMB-BI{IM:3}:IC4_MON')

        
    def state(self):
        
        return "in"
    
    
    def v_position(self):
        
        total = self.v1.value+self.v2.value
        if total>0:
            return (self.v1.value-self.v2.value)/(total)
        else:
            return 0

    def h_position(self):
        
        total = self.h1.value+self.h2.value
        if total>0:
            return (self.h1.value-self.h2.value)/(total)
        else:
            return 0
    
    def reading(self, verbosity=3):
        
        total = self.h1.value + self.h2.value + self.v1.value + self.v2.value
        
        if verbosity>=3:
            print('Reading for {:s} ({:s})'.format(self.name, self.description))
            print('  Horizontal:  {:9.4g}  +  {:9.4g}  =  {:9.4g}'.format(self.h1.value, self.h2.value, self.h1.value+self.h2.value))
            print('    position: {:.3f}'.format(self.h_position()))
            print('  Vertical:    {:9.4g}  +  {:9.4g}  =  {:9.4g}'.format(self.v1.value, self.v2.value, self.v1.value+self.v2.value))
            print('    position: {:.3f}'.format(self.v_position()))

        if verbosity>=2:
            
            print('  Total:  {:9.4g}'.format(total))
            
        return total
    
    
    def current_to_flux(self, current):
        
        energy_keV = self.beam.energy(verbosity=0)
        
        V_ion = 0.036       ## ionization energy of N2 gas in [keV]
        IC_len = 6.0        ## active length of Ion Chamber in [cm]
        qe = 1.602e-19      ## electron charge in [C]

        ## Absorption length [cm] of gas N2 (1 atm, 1.131 g/L) vs E [keV]
        # based on polynomial fit to the calculated abs length data from: henke.lbl.gov/optical_constants/atten2.html 
        # see /home/xf11bm/masa/atten_len_N2* 
        abs_len = 355.21 - 112.26*energy_keV + 11.200*np.square(energy_keV) - 0.10611*np.power(energy_keV,3.0)

        N_abs = current*V_ion/(qe*energy_keV)
        flux = N_abs / (1.0 - np.exp(-IC_len/abs_len))

        return flux
    
    
    def flux(self, verbosity=3):
        
        if self.reading(verbosity=0) < 5e-10:
            return 0.0
        
        h1 = self.current_to_flux(self.h1.value)
        h2 = self.current_to_flux(self.h2.value)
        h_total = h1 + h2
        v1 = self.current_to_flux(self.v1.value)
        v2 = self.current_to_flux(self.v2.value)
        v_total = v1 + v2
        
        total = h_total + v_total
        avg = total*0.5
        
        if verbosity>=3:
            print('Flux for {:s} ({:s})'.format(self.name, self.description))
            print('  Horizontal:  {:9.4g}  +  {:9.4g}  =  {:9.4g} ph/s'.format(h1, h2, h1+h2))
            print('    position: {:.3f}'.format(self.h_position()))
            print('  Vertical:    {:9.4g}  +  {:9.4g}  =  {:9.4g} ph/s'.format(v1, v2, v1+v2))
            print('    position: {:.3f}'.format(self.v_position()))

        if verbosity>=2:
            
            print('  Average:  {:9.4g} ph/s'.format(avg))
            
        return avg 
    
    
    
#ionchamber = IonChamber_CMS(beam=beam)


# CMSBeam
################################################################################
class CMSBeam(object):
    """
    This class represents the 'beam' at the beamline. This collects together aspects
    of querying or changing beam properties, including the energy (or wavelength), 
    the beam intensity (or measuring flux), and so forth.
    """
    
    def __init__(self):
        
        self.mono_bragg_pv = 'XF:11BMA-OP{Mono:DMM-Ax:Bragg}Mtr.RBV'
        
        # (planck constant * speed of light)/(electronic charge)
        self.hc_over_e = 1.23984197e-6 # m^3 kg s^-3 Amp^-1 = eV*m
        self.hc_over_e_keVpA = self.hc_over_e*1e7 # = 12.4 keV*Angstrom
        
        # DMM bilayer pitch in Angstroms, according to Rigaku metrology report
        self.dmm_dsp = 20.1 # Angstroms
        
        
        
        self.mono = BeamlineElement('monochromator', 27.0)
        def transmission(verbosity=0):
            return 1e-7
        self.mono.transmission = transmission

        
        self.attenuator = BeamlineElement('attenuator', 50.0, description="Attenuator/filter box")
        self.attenuator.has_flux = False
        def reading(verbosity=0):
            return self.transmission(verbosity=verbosity)
        self.attenuator.reading = reading
        self.attenuator.transmission = self.transmission

        if False:
            self.fs1 = DiagnosticScreen( 'fs1', 27.5, pv='XF:11BMA-BI{FS:1}', epics_signal=StandardProsilica('XF:11BMA-BI{FS:1-Cam:1}', name='fs1') )
            #self.fs2 = DiagnosticScreen( 'fs2', 29.1, pv='XF:11BMA-BI{FS:2}', epics_signal=StandardProsilica('XF:11BMA-BI{FS:2-Cam:1}', name='fs2') )
            self.fs3 = DiagnosticScreen( 'fs3', 54.0, pv='XF:11BMB-BI{FS:3}', epics_signal=StandardProsilica('XF:11BMB-BI{FS:3-Cam:1}', name='fs3') )
            self.fs4 = DiagnosticScreen( 'fs4', 56.0, pv='XF:11BMB-BI{FS:4}', epics_signal=StandardProsilica('XF:11BMB-BI{FS:4-Cam:1}', name='fs4') )
            self.fs5 = DiagnosticScreen( 'fs5', 70.0, pv='XF:11BMB-BI{FS:Test-Cam:1}', epics_signal=StandardProsilica('XF:11BMB-BI{FS:4-Cam:1}', name='fs5') )
        else:
            # Rely on the fact that these are defined in 20-area-detectors.py
            self.fs1 = DiagnosticScreen( 'fs1', 27.5, pv='XF:11BMA-BI{FS:1}', epics_signal=fs1 )
            #self.fs2 = DiagnosticScreen( 'fs2', 29.1, pv='XF:11BMA-BI{FS:2}', epics_signal=fs2 )
            self.fs3 = DiagnosticScreen( 'fs3', 54.0, pv='XF:11BMB-BI{FS:3}', epics_signal=fs3 )
            self.fs4 = DiagnosticScreen( 'fs4', 56.0, pv='XF:11BMB-BI{FS:4}', epics_signal=fs4 )
            self.fs5 = DiagnosticScreen( 'fs5', 70.0, pv='XF:11BMB-BI{FS:Test-Cam:1}', epics_signal=fs5 )
            
            
        self.bim3 = IonChamber_CMS(beam=self)
        self.beam_defining_slit = s4
        self.bim6 = PointDiode_CMS()
        
        self.GVdsbig = GateValve('GV ds big', 59.0, pv='XF:11BMB-VA{Chm:Det-GV:1}')
        
        
        
        self.elements = []
        
        # Front End
        self.elements.append(ThreePoleWiggler())
        #SR:C03-EPS{PLC:1}Sts:BM_BMPS_Opn-Sts BMPS
        self.elements.append(GateValve('GV1', 20.0, pv='FE:C03A-VA{GV:1}DB:'))
        self.elements.append(GateValve('GV2', 21.0, pv='FE:C03A-VA{GV:2}DB:'))
        
        
        # FOE
        self.elements.append(Shutter('FE shutter', 25.0, pv='XF:11BM-PPS{Sh:FE}'))
        self.elements.append(GateValve('GV', 26.0, pv='FE:C11B-VA{GV:2}'))
        self.elements.append(self.mono)
        self.elements.append(self.fs1)
        # bim1
        # slit0
        # bim2
        self.elements.append(GateValve('GV', 28.0, pv='XF:11BMA-VA{Slt:0-GV:1}'))
        self.elements.append(BeamlineElement('mirror', 29.0))
        self.elements.append(GateValve('GV', 29.0, pv='XF:11BMA-VA{Mir:Tor-GV:1}'))
        self.elements.append(BeamlineElement('fs2 (manual)', 29.1)) # self.elements.append(self.fs2)
        self.elements.append(Shutter('photon shutter', 30.0, pv='XF:11BMA-PPS{PSh}'))
        self.elements.append(GateValve('GV', 30.1, pv='XF:11BMA-VA{PSh:1-GV:1}'))
        
        # Endstation
        self.elements.append(self.bim3)
        # Experimental shutter 49.5
        self.elements.append(self.attenuator)
        self.elements.append(self.fs3)
        self.elements.append(BeamlineElement('KB mirrors', 55.0))
        self.elements.append(self.fs4)
        # im4
        #self.elements.append(GateValve('GV us small', 57.0, pv='XF:11BMB-VA{Slt:4-GV:1}'))
        
        
        self.elements.append(BeamlineElement('sample', 58.0))
        self.elements.append(self.bim6) # dsmon
        self.elements.append(BeamlineElement('WAXS detector', 58.4))
        self.elements.append(self.GVdsbig)
        self.elements.append(BeamlineElement('SAXS detector', 58+5))
        
        
        
        # Sort by position along the beam
        self.elements.sort(key=lambda o: o.zposition, reverse=False)
    
    
    # Monochromator
    ########################################
    
    def energy(self, verbosity=3):
        """
        Returns the current x-ray photon energy (in keV).
        """
        
        # Current angle of monochromator multilayer crystal
        Bragg_deg = caget(self.mono_bragg_pv)
        Bragg_rad = np.radians(Bragg_deg)
        
        wavelength_A = 2.*self.dmm_dsp*np.sin(Bragg_rad)
        wavelength_m = wavelength_A*1e-10

        energy_eV = self.hc_over_e/wavelength_m
        energy_keV = energy_eV/1000.
        
        if verbosity>=3:
            print('E = {:.2f} keV, wavelength = {:.4f} Å, Bragg = {:.6f} rad = {:.4f} deg'.format(energy_keV, wavelength_A, Bragg_rad, Bragg_deg))
            
        elif verbosity>=1:
            print('E = {:.3f} keV'.format(energy_keV))
        
        return energy_keV
        
        
    def wavelength(self, verbosity=3):
        """
        Returns the current x-ray photon wavelength (in Angstroms).
        """
        
        # Current angle of monochromator multilayer crystal
        Bragg_deg = caget(self.mono_bragg_pv)
        Bragg_rad = np.radians(Bragg_deg)
        
        wavelength_A = 2.*self.dmm_dsp*np.sin(Bragg_rad)
        wavelength_m = wavelength_A*1e-10

        # (planck constant * speed of light)/(electronic charge)
        
        energy_eV = self.hc_over_e/wavelength_m
        energy_keV = energy_eV/1000.
        
        if verbosity>=3:
            print('wavelength = {:.4f} Å, E = {:.2f} keV, Bragg = {:.6f} rad = {:.4f} deg'.format(wavelength_A, energy_keV, Bragg_rad, Bragg_deg))
            
        elif verbosity>=1:
            print('wavelength = {:.5f} Å'.format(wavelength_A))
        
        return wavelength_A
    
    
    def setEnergy(self, energy_keV, verbosity=3):
        """
        Set the x-ray beam to the specified energy (by changing the
        monochromator angle.
        """
        
        energy_eV = energy_keV*1000.
        wavelength_m = self.hc_over_e/energy_eV
        wavelength_A = wavelength_m*1.e10
        
        self.setWavelength(wavelength_A, verbosity=verbosity)
        
        return self.energy(verbosity=0)
    
    
    def setWavelength(self, wavelength_A, verbosity=3):
        """
        Set the x-ray beam to the specified wavelength (by changing the
        monochromator angle.
        """
        
        Bragg_deg_initial = caget(self.mono_bragg_pv)
        wavelength_m = wavelength_A*1.e-10
        Bragg_rad = np.arcsin(wavelength_A/(2.*self.dmm_dsp))
        Bragg_deg = np.degrees(Bragg_rad)
        
        print('mono_bragg will move to {:.4f}g deg'.format(Bragg_deg))
        response = input('    Are you sure? (y/[n]) ')
        if response is 'y' or response is 'Y':
            
            mov(mono_bragg, Bragg_deg)
            
            if verbosity>=1:
                print('mono_bragg moved from {.4f} deg to {.4f} deg'.format(Bragg_deg_initial, Bragg_deg))
        
        elif verbosity>=1:
            print('No move was made.')
            
        return self.wavelength(verbosity=verbosity)

    
    # Slits
    ########################################
    
    def size(self, verbosity=3):
        """
        Returns the current beam size (rough estimate).
        The return is (size_horizontal, size_vertical) (in mm).
        """
        size_h = self.beam_defining_slit.xg.user_readback.value
        size_v = self.beam_defining_slit.yg.user_readback.value
        
        if verbosity>=3:
            print('Beam size:')
            print('  horizontal = {:.3f} mm'.format(size_h))
            print('  vertical   = {:.3f} mm'.format(size_v))
        
        return size_h, size_v

    
    def setSize(self, horizontal, vertical, verbosity=3):
        """
        Sets the beam size.
        """
        
        h, v = self.size(verbosity=0)
        
        if verbosity>=3:
            print('Changing horizontal beam size from {:.3f} mm to {:.3f} mm'.format(h, horizontal))
        self.beam_defining_slit.xg.user_setpoint.value = horizontal
        
        if verbosity>=3:
            print('Changing vertical beam size from {:.3f} mm to {:.3f} mm'.format(v, vertical))
        
        self.beam_defining_slit.yg.user_setpoint.value = vertical
    
    
    def divergence(self, verbosity=3):
        """
        Returns the beamline divergence.
        This is based on the Front End (FE) slits. The return is
        (horizontal, vertical) (in mrad).
        """
        
        distance_m = 10.0 # distance from source to slits
        
        horizontal_mm = caget('FE:C11B-OP{Slt:12-Ax:X}t2.C')
        vertical_mm = caget('FE:C11B-OP{Slt:12-Ax:Y}t2.C')
        
        horizontal_mrad = horizontal_mm/distance_m
        vertical_mrad = vertical_mm/distance_m
        
        if verbosity>=3:
            print('Beam divergence:')
            print('  horizontal = {:.3f} mrad'.format(horizontal_mrad))
            print('  vertical   = {:.3f} mrad'.format(vertical_mrad))
        
        return horizontal_mrad, vertical_mrad
        
    
    def setDivergence(self, horizontal, vertical, verbosity=3):
        """
        Set beamline divergence (in mrad).
        This is controlled using the Front End (FE) slits.
        """
        
        h, v = self.divergence(verbosity=0)

        distance_m = 10.0 # distance from source to slits
        
        horizontal_mm = horizontal*distance_m
        vertical_mm = vertical*distance_m
        
        if horizontal<0:
            if verbosity>=1:
                print("Horizontal divergence less than zero ({}) doesn't make sense.".format(horizontal))
            
        elif horizontal>1.5:
            if verbosity>=1:
                print("Horizontal divergence should be less than 1.5 mrad.")
                
        else:
            if verbosity>=3:
                print('Changing horizontal divergence from {:.3f} mrad to {:.3f} mrad.'.format(h, horizontal))
            caput('FE:C11B-OP{Slt:12-Ax:X}size', horizontal_mm)
        
        
        if vertical<0:
            if verbosity>=1:
                print("Vertical divergence less than zero ({}) doesn't make sense.".format(vertical))
            
        elif vertical>0.15:
            if verbosity>=1:
                print("Vertical divergence should be less than 0.15 mrad.")
                
        else:
            if verbosity>=3:
                print('Changing vertical divergence from {:.3f} mrad to {:.3f} mrad.'.format(v, vertical))
            caput('FE:C11B-OP{Slt:12-Ax:Y}size', vertical_mm)
        

    
    # Experimental Shutter
    ########################################
    
    def is_on(self, verbosity=3):
        '''Returns true if the beam is on (experimental shutter open).'''
        
        blade1 = caget('XF:11BMB-OP{PSh:2}Pos:1-Sts')
        blade2 = caget('XF:11BMB-OP{PSh:2}Pos:2-Sts')
        
        if blade1==1 and blade2==1:
            if verbosity>=4:
                print('Beam on (shutter open).')
            
            return True
        
        else:
            if verbosity>=4:
                print('Beam off (shutter closed).')
            
            return False
    
    
    def on(self, verbosity=3, wait_time=0.005, poling_period=0.10, retry_time=2.0, max_retries=5):
        '''Turn on the beam (open experimental shutter).'''
        
        if self.is_on(verbosity=0):
            if verbosity>=4:
                print('Beam on (shutter already open.)')
                
        else:
            
            itry = 0
            while (not self.is_on(verbosity=0)) and itry<max_retries:
            
                # Trigger the shutter to toggle state
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M112=1')
                sleep(wait_time)
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M111=1')
                sleep(wait_time)
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M112=0')
                sleep(wait_time)
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M111=1')
                sleep(wait_time)
                
                # Give the system a chance to update
                start_time = time.time()
                while (not self.is_on(verbosity=0)) and (time.time()-start_time)<retry_time:
                    if verbosity>=5:
                        print('  try {:d}, t = {:02.2f} s, state = {:s}'.format(itry+1, (time.time()-start_time), 'OPEN_____' if self.is_on(verbosity=0) else 'CLOSE===='))
                    sleep(poling_period)
                
                itry += 1
                

            if verbosity>=4:
                if self.is_on(verbosity=0):
                    print('Beam on (shutter opened).')
                else:
                    print("Beam off (shutter didn't open).")

    
    def off(self, verbosity=3, wait_time=0.005, poling_period=0.10, retry_time=2.0, max_retries=5):
        '''Turn off the beam (close experimental shutter).'''
        
        if self.is_on(verbosity=0):
            
            itry = 0
            while self.is_on(verbosity=0) and itry<max_retries:
                # Trigger the shutter to toggle state
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M112=1')
                sleep(wait_time)
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M111=1')
                sleep(wait_time)
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M112=0')
                sleep(wait_time)
                caput('XF:11BMB-CT{MC:06}Asyn.AOUT','M111=1')
                sleep(wait_time)

                # Give the system a chance to update
                start_time = time.time()
                while self.is_on(verbosity=0) and (time.time()-start_time)<retry_time:
                    if verbosity>=5:
                        print('  try {:d}, t = {:02.2f} s, state = {:s}'.format(itry+1, (time.time()-start_time), 'OPEN_____' if self.is_on(verbosity=0) else 'CLOSE===='))
                    sleep(poling_period)
                
                itry += 1



            if verbosity>=4:
                if self.is_on(verbosity=0):
                    print("Beam on (shutter didn't close).")
                else:
                    print('Beam off (shutter closed).')
                
        else:
            if verbosity>=4:
                print('Beam off (shutter already closed).')
            
    
    
    # Attenuator/Filter Box
    ########################################

    def transmission(self, verbosity=3):
        """
        Returns the current beam transmission through the attenuator/filter box.
        To change the transmission, use 'setTransmission'.
        """
        
        energy_keV = self.energy(verbosity=0)
        
        if energy_keV < 6.0 or energy_keV > 18.0:
            print('Transmission data not available at the current X-ray energy ({.2f} keV).'.format(energy_keV))
            
        else:
            
            # The states of the foils in the filter box
            N = [ caget('XF:11BMB-OP{{Fltr:{:d}}}Pos-Sts'.format(ifoil)) for ifoil in range(1, 8+1) ]
            
            tr_tot = self.calc_transmission_filters(N, verbosity=verbosity)
                    
            return tr_tot


    def calc_transmission_filters(self, filter_settings, energy_keV=None, verbosity=3):
        """
        Returns the computed transmission value for the given configuration of
        foils. Note that the foils are not actually moved. This is just a
        calculation.
        
        Parameters
        ----------
        filter_settings : array of length 8
            Each element must be either a zero (foil removed) or a 1 (foil blocking 
            beam)
        energy_keV : float
            If 'None', the current energy is used. If specified, the calculation 
            is performed for the requested energy.
            
        Returns                
                

        -------
        transmission : float
            The computed transmission value of the x-ray beam through the filter box.
        """
        
        if energy_keV is None:
            energy_keV = self.energy(verbosity=0)
            
        if len(filter_settings) != 8:
            print('States for all eight foils must be specified.')

        else:
            N = filter_settings
            
            E = energy_keV
            E2 = np.square(E)
            E3 = np.power(E, 3)
                        

            # foil thickness blocking the beam
            N_Al = N[0] + 2*N[1] + 4*N[2] + 8*N[3]
            N_Nb = N[4] + 2*N[5] + 4*N[6] + 8*N[7]

            d_Nb = 0.1      # Thickness [mm] of one Nb foil 
            d_Al = 0.25     # Thickness [mm] of one Al foil 

            # Absorption length [mm] based on fits to LBL CXRO data for 6 < E < 19 keV
            l_Nb = 1.4476e-3 - 5.6011e-4 * E + 1.0401e-4 * E2 + 8.7961e-6 * E3
            l_Al = 5.2293e-3 - 1.3491e-3 * E + 1.7833e-4 * E2 + 1.4001e-4 * E3

            # transmission factors
            tr_Nb = np.exp(-N_Nb*d_Nb/l_Nb)
            tr_Al = np.exp(-N_Al*d_Al/l_Al)
            tr_tot = tr_Nb*tr_Al
                
            if verbosity>=5:
                print('  state:      {} T = {:.6g}'.format(filter_settings, tr_tot))
            if verbosity>=4:
                print('{:d} × 0.25 mm Al ({:.4g}) and {:d} × 0.10 mm Nb ({:.4g})'.format(N_Al, tr_Al, N_Nb, tr_Nb) )
            if verbosity>=1:
                print('transmission = {:.6g}'.format(tr_tot))
                
            return tr_tot
            
        

    def set_attenuation_filters(self, filter_settings, verbosity=3):
        """
        Sets the positions (in/out) for each of the foils in the attenuator/
        filter box. The input 'filter_settings' should be an array of length
        8, where each element is either a zero (foil removed) or a 1 (foil
        blocking beam).
        """
        
        if verbosity>=4:
            print('Filters:')
            # The states of the foils in the filter box
            filters_initial = [ caget('XF:11BMB-OP{{Fltr:{:d}}}Pos-Sts'.format(ifoil)) for ifoil in range(1, 8+1) ]
            print('  initial:    {} T = {:.6g}'.format(filters_initial, self.calc_transmission_filters(filters_initial, verbosity=0)))
            print('  requested:  {} T = {:.6g}'.format(filter_settings, self.calc_transmission_filters(filter_settings, verbosity=0)))
        
        if len(filter_settings) != 8:
            print('States for all eight foils must be specified.')
            
        else:
            
            for i, state in enumerate(filter_settings):
                
                ifoil = i+1
                
                if state==1:
                    # Put foil #ifoil into the beam
                    caput( 'XF:11BMB-OP{{Fltr:{:d}}}Cmd:In-Cmd'.format(ifoil) , 1 )
                    
                elif state==0:
                    # Remove foil #ifoil
                    caput( 'XF:11BMB-OP{{Fltr:{:d}}}Cmd:Out-Cmd'.format(ifoil) , 1 )
                    
                else:
                    if verbosity>=3:
                        state_actual = caget( 'XF:11BMB-OP{{Fltr:{:d}}}Pos-Sts'.format(ifoil) )
                        state_actual_str = 'IN' if state_actual is 1 else 'OUT'
                        print('WARNING: Filter state {} not recognized. Filter {:d} is {:s}.'.format(state, ifoil, state_actual_str))
                    

        
            sleep(1.) # Wait for filter box to settle
            
        if verbosity>=4:
            filters_final = [ caget('XF:11BMB-OP{{Fltr:{:d}}}Pos-Sts'.format(ifoil)) for ifoil in range(1, 8+1) ]
            print('  final:      {} T = {:.6g}'.format(filters_final, self.calc_transmission_filters(filters_final, verbosity=0)))

        
    def setTransmission(self, transmission, verbosity=3):
        """
        Sets the transmission through the attenuator/filter box.
        Because the filter box has a discrete set of foils, it is impossible to
        exactly match a given transmission value. A nearby value will be
        selected.
        """
        
        energy_keV = self.energy(verbosity=0)
        
        if energy_keV < 6.0 or energy_keV > 18.0:
            print('Transmission data not available at the current X-ray energy ({.2f} keV).'.format(energy_keV))
            
        elif transmission > 1.0:
            print('A transmission above 1.0 is not possible.')
            
        elif transmission < 1e-10:
            print('A transmission this low ({:g}) cannot be reliably achieved.'.format(transmission))
            
        else:
            
            E = energy_keV
            E2 = np.square(E)
            E3 = np.power(E, 3)
            
            d_Nb = 0.1      # Thickness [mm] of one Nb foil 
            d_Al = 0.25     # Thickness [mm] of one Al foil 

            # Absorption length [mm] based on fits to LBL CXRO data for 6 < E < 19 keV
            l_Nb = 1.4476e-3 - 5.6011e-4 * E + 1.0401e-4 * E2 + 8.7961e-6 * E3
            l_Al = 5.2293e-3 - 1.3491e-3 * E + 1.7833e-4 * E2 + 1.4001e-4 * E3

            d_l_Nb = d_Nb/l_Nb
            d_l_Al = d_Al/l_Al

            # Number of foils to be inserted (equivalent to "XIA_attn.mac" from X9) 
            #N_Nb = int(-log(transmission)/d_l_Nb)
            ##N_Al = int((-log(transmission) - N_Nb*d_l_Nb)/(d_l_Al-0.5))
            #N_Al = int((-log(transmission) - N_Nb*d_l_Nb)/d_l_Al)

            # Number of foils to be inserted (picks a set that gives smallest deviation from requested transmission)
            dev = []
            for i in np.arange(16):
                for j in np.arange(16):
                    dev_ij = abs(transmission - exp(-i*d_l_Nb)*exp(-j*d_l_Al))
                    dev.append(dev_ij)
                    if (dev_ij == min(dev)):
                        N_Nb = i                    # number of Nb foils selected
                        N_Al = j                    # number of Al foils selected
                
                


            N = []
            state = N_Al
            for i in np.arange(4):
                N.append(state % 2)
                state = int(state/2)

            state = N_Nb
            for i in np.arange(4):
                N.append(state % 2)
                state = int(state/2)

            self.set_attenuation_filters(N, verbosity=verbosity)

        
        return self.transmission(verbosity=verbosity)





    # Flux estimates at various points along the beam
    ########################################
    
    # TBD
    
    
    # Flux diagnostics
    ########################################
    
    def fluxes(self, verbosity=3):
        """
        Outputs a list of fluxes at various points along the beam. Also checks 
        the state (in or out of the beam) of various components, to help identify
        if anything could be blocking the beam.
        """
        
        if verbosity>=1:
            print('+--------+------------------+-----+-------------+-------------+-------------+')
            print('| pos    | name             |path | reading     | flux (ph/s) | expected    |')
            print('|--------|------------------|-----|-------------|-------------|-------------|')
            
        
        last_z = -100
        beam = True
        
        flux_expected = None
        
        for element in self.elements:
            
            state = element.state()
            if state is 'block':
                beam = False
            
            if verbosity>=4:
                if element.zposition >= 0 and last_z < 0:
                    print('| Front End                 |     |             |             |             |')
                if element.zposition > 25 and last_z < 25:
                    print('| FOE                       |     |             |             |             |')
                if element.zposition > 50 and last_z < 50:
                    print('| Endstation                |     |             |             |             |')
            last_z = element.zposition
            flux_expected
            if verbosity>=1:
                
                
                if state is 'in':
                    if beam:
                        path = '(|)'
                    else:
                        path = '(-)'
                elif state is 'out':                
                

                    if beam:
                        path = ' | '
                    else:
                        path = '---'
                elif state is 'block':
                    path = '[X]'
                    beam = False
                
                elif state is 'undefined':
                    if beam:
                        path = '?|?'
                    else:
                        path = '?-?'
                
                else:
                    path = '???'
                    


                
                
                if flux_expected is None or not beam:
                    flux_expected_str = ''
                else:
                    flux_expected_str = '{:11.3g}'.format(flux_expected)
                    flux_expected *= element.transmission(verbosity=0)


                
                
                if callable(getattr(element, 'reading', None)):
                    reading_str = '{:11.3g}'.format(element.reading(verbosity=0))
                    state = element.state()
                    if element.has_flux and (state=='in' or state=='block'):
                        flux_cur = element.flux(verbosity=0)
                        flux_expected = flux_cur
                        flux_str = '{:11.3g}'.format(flux_cur)
                    else:
                        flux_str = ''
                    
                else:
                    reading_str = ''
                    flux_str = ''
                    
            
                
                print('|{:5.1f} m | {:16.16} | {:s} | {:11.11} | {:11.11} | {:11.11} |'.format(element.zposition, element.name, path, reading_str, flux_str, flux_expected_str))
                
                
            #beam = True # For testing
                
                
        if verbosity>=1:
            print('+--------+------------------+-----+-------------+-------------+-------------+')

            
            


    # End class CMSBeam(object)
    ########################################
    


beam = CMSBeam()


class Beamline(object):
    '''Generic class that encapsulates different aspects of the beamline.
    The intention for this object is to have methods that activate various 'standard'
    protocols or sequences of actions.'''

    def __init__(self, **kwargs):
        
        self.md = {}
        self.current_mode = 'undefined'
        
        
    def mode(self, new_mode):
        '''Tells the instrument to switch into the requested mode. This may involve
        moving detectors, moving the sample, enabling/disabling detectors, and so
        on.'''
        
        getattr(self, 'mode'+new_mode)()
        
        
    def get_md(self, prefix=None, **md):
        '''Returns a dictionary of the current metadata.
        The 'prefix' argument is prepended to all the md keys, which allows the
        metadata to be grouped with other metadata in a clear way. (Especially,
        to make it explicit that this metadata came from the beamline.)'''
        
        # Update internal md
        #self.md['key'] = value

        md_return = self.md.copy()
    
        # Add md that may change
        md_return['mode'] = self.current_mode
    
        # Include the user-specified metadata
        md_return.update(md)

        # Add an optional prefix
        if prefix is not None:
            md_return = { '{:s}{:s}'.format(prefix, key) : value for key, value in md_return.items() }
    
        return md_return
            
        
    def comment(self, text, logbooks=None, tags=None, append_md=True, **md):
        
        text += '\n\n[comment for beamline: {}]'.format(self.__class__.__name__)
        
        if append_md:
        
            # Global md
            md_current = { k : v for k, v in RE.md.items() }
            
            # Beamline md
            md_current.update(self.get_md())
            
            # Specified md
            md_current.update(md)
            
            text += '\n\n\nMetadata\n----------------------------------------'
            for key, value in sorted(md_current.items()):
                text += '\n{}: {}'.format(key, value)
        
        logbook.log(text, logbooks=logbooks, tags=tags)
        
        
    def log_motors(self, motors, verbosity=3, **md):
      
        log_text = 'Motors\n----------------------------------------\nname | position | offset | direction |\n'
      
        for motor in motors:
            offset = float(caget(motor.prefix+'.OFF'))
            direction = int(caget(motor.prefix+'.DIR'))
            log_text += '{} | {} | {} | {} |\n'.format(motor.name, motor.user_readback.value, offset, direction)
      
      
        md_current = { k : v for k, v in RE.md.items() }
        md_current.update(md)
        log_text += '\nMetadata\n----------------------------------------\n'
        for k, v in sorted(md_current.items()):
            log_text += '{}: {}\n'.format(k, v)
            
        if verbosity>=3:
            print(log_text)
            
        self.comment(log_text)
            

                



class CMS_Beamline(Beamline):
    '''This object collects together various standard protocols and sequences
    of action used on the CMS (11-BM) beamline at NSLS-II.'''
    
    
    def __init__(self, **kwargs):
        
        super().__init__(**kwargs)
        
        self.beam = beam
        
        from epics import PV
        
        self._chamber_pressure_pv = PV('XF:11BMB-VA{Chm:Det-TCG:1}P-I')
        
    
    
    def modeAlignment(self, verbosity=3):
        
        self.current_mode = 'undefined'
        
        # TODO: Check what mode (TSAXS, GISAXS) and respond accordingly
        # TODO: Check if gate valves are open and flux is okay (warn user)
        
        
        self.beam.off()
        self.beam.setTransmission(5e-4)
        
        #mov( [DETx, DETy], [0, 0] )
        self.beam.bim6.insert()
        
        caput('XF:11BMB-BI{IM:2}EM180:Acquire', 1) # Turn on bim6
        
        self.current_mode = 'alignment'
        
        self.beam.bim6.reading()
        
        
        
    def modeMeasurement(self, verbosity=3):
        
        self.current_mode = 'undefined'
        
        self.beam.off()
        self.beam.setTransmission(1)
        
        #mov(DETy, -16)
        self.beam.bim6.retract()
        
        caput('XF:11BMB-BI{IM:2}EM180:Acquire', 0) # Turn off bim6
        
        self.current_mode = 'measurement'
        
        # Check if gate valves are open
        if self.beam.GVdsbig.state() is not 'out' and verbosity>=1:
            print('Warning: Sample chamber gate valve (large, downstream) is not open.')
            
        
        
        
    def modeBeamstopAlignment(self, verbosity=3):
        '''Places bim6 (dsmon) as a temporary beamstop.'''
        
        mov(DETy, -6.1)
        
        
        
    def ventChamber(self, verbosity=3):
        
        
        # Close large gate valve (downstream side of sample chamber)
        caput('XF:11BMB-VA{Chm:Det-GV:1}Cmd:Cls-Cmd',1)

        # Close small gate valve (upstream side of sample chamber)
        #caput('XF:11BMB-VA{Slt:4-GV:1}Cmd:Cls-Cmd',1)

        # Close valve connecting sample chamber to vacuum pump
        caput('XF:11BMB-VA{Chm:Det-IV:1}Cmd:Cls-Cmd',1)
        
        sleep(0.5)
        
        # Soft-open the upstream vent-valve
        caput('XF:11BMB-VA{Chm:Smpl-VV:1}Cmd:Cls-Cmd', 1)
        sleep(1.0)
        caput('XF:11BMB-VA{Chm:Smpl-VV:1_Soft}Cmd:Opn-Cmd', 1)
        
        
        
        self.chamberPressure(range_high=100)
        
        # Fully open the upstream vent-vale
        caput('XF:11BMB-VA{Chm:Smpl-VV:1_Soft}Cmd:Cls-Cmd', 1)
        sleep(1.0)
        caput('XF:11BMB-VA{Chm:Smpl-VV:1}Cmd:Opn-Cmd', 1)

        # Fully open the downstream vent-vale
        caput('XF:11BMB-VA{Chm:Det-VV:1_Soft}Cmd:Cls-Cmd', 1)
        sleep(1.0)
        caput('XF:11BMB-VA{Chm:Det-VV:1}Cmd:Opn-Cmd', 1)
        
        self.chamberPressure(range_high=1000)
        
        if verbosity>=1:
            print('Sample chamber is ready to be opened.')
        
        
        
    def chamberPressure(self, range_low=None, range_high=None, readout_period=1.0, verbosity=3):
        '''Monitors the pressure in the sample chamber, printing the current value.
        If range arguments are provided, the monitoring will end once the pressure
        is outside the range.
        '''
        
        monitor = True
        while monitor:
            
            try:
            
                if range_low is not None and self._chamber_pressure_pv.get()<range_low:
                    monitor = False
                    
                if range_high is not None and self._chamber_pressure_pv.get()>range_high:
                    monitor = False

                P_mbar = self._chamber_pressure_pv.get()
                P_atm = P_mbar*0.000986923
                P_torr = P_mbar*0.750062
                P_kPa = P_mbar*0.1
                P_psi = P_mbar = 0.0145038

                if verbosity>=4:
                    print('Sample chamber pressure: {:8.2f} mbar = {:5.3f} atm = {:7.3f} torr = {:4.1g} kPa     \r'.format(P_mbar, P_atm, P_torr, P_kPa), end='', flush=True)
                elif verbosity>=2:
                    print('Sample chamber pressure: {:8.2f} mbar ({:5.3f} atm)    \r'.format(P_mbar, P_atm), end='', flush=True)
                    
                sleep(readout_period)
                
                
            except KeyboardInterrupt:
                monitor = False
                
        
        
    def pumpChamber(self, readout_delay=0.2):
        
        
        # Close vent-valves
        caput('XF:11BMB-VA{Chm:Smpl-VV:1_Soft}Cmd:Cls-Cmd', 1)
        sleep(0.5)
        caput('XF:11BMB-VA{Chm:Smpl-VV:1}Cmd:Cls-Cmd', 1)
        sleep(0.5)
        caput('XF:11BMB-VA{Chm:Det-VV:1_Soft}Cmd:Cls-Cmd', 1)
        sleep(0.5)
        caput('XF:11BMB-VA{Chm:Det-VV:1}Cmd:Cls-Cmd', 1)
        sleep(0.2)
        
        # Turn on pump (if necessary)
        if caget('XF:11BMB-VA{Chm:Det-Pmp:1}Sts:Enbl-Sts')==0:
            caput('XF:11BMB-VA{Chm:Det-Pmp:1}Cmd:Enbl-Cmd', 0)
            sleep(0.2)
            caput('XF:11BMB-VA{Chm:Det-Pmp:1}Cmd:Enbl-Cmd', 1)
        
        # Soft-open valve to pump
        caput('XF:11BMB-VA{Chm:Det-IV:1}Cmd:Cls-Cmd', 1)
        sleep(1.0)
        caput('XF:11BMB-VA{Chm:Det-IV:1_Soft}Cmd:Opn-Cmd', 1)
        sleep(0.2)
        
        sleep(5.0)
        # Check pump again
        if caget('XF:11BMB-VA{Chm:Det-Pmp:1}Sts:Enbl-Sts')==0:
            caput('XF:11BMB-VA{Chm:Det-Pmp:1}Cmd:Enbl-Cmd', 0)
            sleep(0.2)
            caput('XF:11BMB-VA{Chm:Det-Pmp:1}Cmd:Enbl-Cmd', 1)
        
        
        self.chamberPressure(range_low=500)

        # Fully open valve to pump
        caput('XF:11BMB-VA{Chm:Det-IV:1_Soft}Cmd:Cls-Cmd', 1)
        sleep(1.0)
        caput('XF:11BMB-VA{Chm:Det-IV:1}Cmd:Opn-Cmd', 1)
        sleep(0.2)
        
        self.chamberPressure(range_low=200)
        
        
    def openChamberGateValve(self):
        
        caput('XF:11BMB-VA{Chm:Det-GV:1}Cmd:Opn-Cmd', 1) # Large (downstream)
        #caput('XF:11BMB-VA{Slt:4-GV:1}Cmd:Opn-Cmd',1) # Small (upstream)


    def closeChamberGateValve(self):
        
        caput('XF:11BMB-VA{Chm:Det-GV:1}Cmd:Cls-Cmd', 1) # Large (downstream)
        #caput('XF:11BMB-VA{Slt:4-GV:1}Cmd:Cls-Cmd',1) # Small (upstream)
        
        
    # Metatdata methods
    ########################################
    
    def get_md(self, prefix=None, **md):
        
        md_current = self.md.copy()
        md_current['calibration_energy_keV'] = self.beam.energy(verbosity=0)
        md_current['calibration_wavelength_A'] = self.beam.wavelength(verbosity=0)
        
        h, v = self.beam.size(verbosity=0)
        md_current['beam_size_x_mm'] = h
        md_current['beam_size_y_mm'] = v
        h, v = self.beam.divergence(verbosity=0)
        md_current['beam_divergence_x_mrad'] = h
        md_current['beam_divergence_y_mrad'] = v
        
        md_current['beamline_mode'] = self.current_mode
        
        md_current.update(md)
        
        # Add an optional prefix
        if prefix is not None:
            md_current = { '{:s}{:s}'.format(prefix, key) : value for key, value in md_current.items() }
        
        return md_current
    

    def setMetadata(self, verbosity=3):
        '''Guides the user through setting some of the required and recommended
        meta-data fields.'''
        
        if verbosity>=3:
            print('This will guide you through adding some meta-data for the upcoming experiment.')
        if verbosity>=4:
            print('You can accept default values (shown in square [] brackets) by pressing enter. You can leave a value blank (or put a space) to skip that entry.')


        # Set some values automatically
        month = int(time.strftime('%m'))
        if month<=4:
            cycle = 1
        elif month<=8:
            cycle = 2
        else:
            cycle = 3    
        RE.md['experiment_cycle'] = '{:s}-{:d}'.format( time.strftime('%Y'), cycle )
        
        RE.md['calibration_energy_keV'] = self.beam.energy(verbosity=0)
        RE.md['calibration_wavelength_A'] = self.beam.wavelength(verbosity=0)
        
        # TODO:
        # RE.md['calibration_detector_distance_m'] =
        # RE.md['calibration_detector_x0'] =
        # RE.md['calibration_detector_y0'] = 
        
            
        
        # Ask the user some questions
        
        questions = [
            ['experiment_proposal_number', 'Proposal number'] ,
            ['experiment_SAF_number', 'SAF number'] ,
            ['experiment_group', 'User group (e.g. PI)'] ,
            ['experiment_user', 'The specific user/person running the experiment'] ,
            ['experiment_project', 'Project name/code'] ,
            ['experiment_type', 'Type of experiments/measurements (SAXS, GIWAXS, etc.)'] ,
            ]
        
        
        # TBD:
        # Path where data will be stored?

        self._dialog_total_questions = len(questions)
        self._dialog_question_number = 1
        
        for key, text in questions:
            try:
                self._ask_question(key, text)
            except KeyboardInterrupt:
                return
            
        if verbosity>=4:
            print('You can also add/edit metadata directly using the RE.md object.')
        
        

    def _ask_question(self, key, text, default=None):

        if default is None and key in RE.md:
                default = RE.md[key]
        
        if default is None:
            ret = input('  Q{:d}/{:d}. {:s}: '.format(self._dialog_question_number, self._dialog_total_questions, text) )
        
        else:
            ret = input('  Q{:d}/{:d}. {:s} [{}]: '.format(self._dialog_question_number, self._dialog_total_questions, text, default) )
            if ret=='':
                ret = default
                
            
        if ret!='' and ret!=' ':
            RE.md[key] = ret
            
        self._dialog_question_number += 1
            
        
    # Logging methods
    ########################################
        
    def logAllMotors(self, verbosity=3, **md):
        log_pos()
        
        motor_list = [
                    mono_bragg ,
                    mono_pitch2 ,
                    mono_roll2 ,
                    mono_perp2 ,
                    mir_usx ,
                    mir_dsx ,
                    mir_usy ,
                    mir_dsyi ,
                    mir_dsyo ,
                    mir_bend ,
                    s0.tp ,
                    s0.bt ,
                    s0.ob ,
                    s0.ib ,
                    s1.xc ,
                    s1.xg ,
                    s1.yc ,
                    s1.yg ,
                    s2.xc ,
                    s2.xg ,
                    s2.yc ,
                    s2.yg ,
                    s3.xc ,
                    s3.xg ,
                    s3.yc ,
                    s3.yg ,
                    s4.xc ,
                    s4.xg ,
                    s4.yc ,
                    s4.yg ,
                    s5.xc ,
                    s5.xg ,
                    s5.yc ,
                    s5.yg ,
                    bim3y ,
                    fs3y ,
                    bim4y ,
                    bim5y ,
                    smx ,
                    smy ,
                    sth ,
                    schi ,
                    sphi ,
                    srot ,
                    strans ,
                    camx ,
                    camy ,
                    cam2x ,
                    cam2z ,
                    DETx ,
                    DETy ,
                    WAXSx ,
                    SAXSx ,
                    SAXSy ,
                    bsx , 
                    bsy ,
                    bsphi ,
                    armz ,
                    armx ,
                    armphi ,
                    army ,
                    armr ,
                    ]
        
        self.log_motors(motor_list, verbosity=verbosity, **md)
        
       
    # End class CMSBeam(object)
    ########################################
        
        

cms = CMS_Beamline()

def get_beamline():
    return cms


