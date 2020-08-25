import time
from ophyd import Device
# Moxa ioLogik and MFC controller
# Moxa: AO E1241, IP: 10.11.130.111    Chan1-4
# Moxa: DIO E1214, IP: 10.11.130.112   Chan1-6   with relay
# Moxa: RTD  E1260, IP: 10.11.130.113  Chan1-6   reading only
# Moxa: Thermal Couple  E1241, IP: 10.11.130.114   Chan1-8  reading only
# Moxa: AO E1241, IP: 10.11.130.115   Chan5-8
# Moxa: AI E1240, IP: 10.11.130.116   Chan1-8    reading only 

# MFC: 

# PV list of Moxa ioLogik:: AO, Analog Output
class AOpv(object):    
    def __init__(self, ii):
        self.name = 'AO_Chan{}'.format(ii)
        self.sp = 'XF:11BMB-ES{}AO:{}-SP'.format('{IO}', ii)
        self.sts = 'XF:11BMB-ES{}AO:{}-RB'.format('{IO}', ii)
AO=[None]
for ii in range(1, 9):
    AO.append(AOpv(ii))

# PV list of Moxa ioLogik:: AI, Analog Input
class AIpv(object):    
    def __init__(self, ii):
        self.name = 'AI_Chan{}'.format(ii)
        self.sts = 'XF:11BMB-ES{}AI:{}-I'.format('{IO}', ii)
AI=[None]
for ii in range(1, 9):
    AI.append(AIpv(ii))

# PV list of Moxa ioLogik:: DO(Relay), Digital Output
class Relaypv(object):    
    def __init__(self, ii):
        self.name = 'Relay_Chan{}'.format(ii)
        self.sp = 'XF:11BMB-ES{}DO:{}-Cmd'.format('{IO}', ii)
        self.sts = 'XF:11BMB-ES{}DO:{}-Sts'.format('{IO}', ii)
Relay=[None]
for ii in range(1, 7):
    Relay.append(Relaypv(ii))

# PV list of Moxa ioLogik:: DI, Digital Input
class DIpv(object):    
    def __init__(self, ii):
        self.name = 'DI_Chan{}'.format(ii)
        self.sts = 'XF:11BMB-ES{}DI:{}-Sts'.format('{IO}', ii)
DI=[None]
for ii in range(1, 7):
    DI.append(DIpv(ii))

#PV list of Moxa ioLogik:: RTD
class RTDpv(object):    
    def __init__(self, ii):
        self.name = 'RTD_Chan{}'.format(ii)
        self.sts = 'XF:11BMB-ES{}T:{}-I'.format('{IO:RTD}', ii)
        
RTD=[None]
for ii in range(1, 7):
    RTD.append(RTDpv(ii))

# PV list of Moxa ioLogik:: TC, Thermal Couple
class TCpv(object):    
    def __init__(self, ii):
        self.name = 'TC_Chan{}'.format(ii)
        self.sts = 'XF:11BMB-ES{}T:{}-I'.format('{IO:TC}', ii)
TC=[None]
for ii in range(1, 7):
    TC.append(TCpv(ii))


#Format: ioL.read(RTD[2]), ioL.set(AO[1], 5), ioL.setOn(Relay[1]), ioL.setOff([Relay[2])

class ioLogik(Device):
    
    def __init__(self, prefix='', *args, read_attrs=None, configuration_attrs=None,name='ioLogik', parent=None, **kwargs):

        super().__init__(prefix=prefix, *args, read_attrs=read_attrs, configuration_attrs=configuration_attrs, name=name, parent=parent, **kwargs)
    
    def read(self, port):
        if port is not None and port in AO+AI+RTD+TC+Relay+DI:
            return caget(port.sts)
        else:
            print('The port is not valid')
            
    def set(self, port, val, verbosity=3):
        if port is not None and port in AO:
            if val >10 or val <0:
                print('Out of input range. It must be in range of (0, 10).')
                return caput(port.sp, 0)
            else:
                caput(port.sp, val)
                if verbosity >=3:
                    time.sleep(0.2)
                    return print('The {} is set as {}.'.format(port.name, self.read(port)))
        else:
            print('The port is not valid.')
    
    def setOn(self, port):
        if port is not None and port in Relay:
            caput(port.sp, 1)
            #re-check the port value
            if self.read(port)==1:
                print('{} is turned on.'.format(port.name))
            else:
                print('{} is turned off.'.format(port.name))
        else:
            print('The port is not valid')

    def setOff(self, port):
        if port is not None and port in Relay:
            caput(port.sp, 0)
            time.sleep(.2)
            if self.read(port)==1:
                print('{} is turned on.'.format(port.name))
            else:
                print('{} is turned off.'.format(port.name))
        else:    
            print('The port is not valid')
            
    def readRH(self, AI_chan, temperature=25.0, voltage_supply=5.0, coeff_slope=0.030, coeff_offset=0.787, verbosity=3):
        voltage_out = self.read(AI[AI_chan])
        corr_voltage_out = voltage_out * (5.0 / voltage_supply)
        #For sensor #220 used for SVA chamber
        #coeff_offset = 0.788 #from the certificate
        #coeff_offset = 0.746 #from the environment of RH=0
        #coeff_slope = 0.029

        #For sensor used for Linkam tensile stage
        #coeff_offset = 0.787
        #coeff_slope = 0.030

        #For sensor 114 used for environmental bar
        #coeff_offset = 0.787
        #coeff_slope = 0.030

        #For sensor 43 used in humidity stage
        coeff_offset = 0.816887
        coeff_slope = 0.028813

        sensor_RH = (corr_voltage_out - coeff_offset) / coeff_slope

        true_RH = sensor_RH / (1.0546 - 0.00216 * temperature)      # T in [degC]
        
        if verbosity >= 3:
            print('Raw sensor RH = {:.3f} pct.'.format(sensor_RH))
            print('T-corrected RH = {:.3f} pct at {:.3f} degC.'.format(true_RH, temperature))
        return true_RH


# MFC PVs: 
#FlowRate_Sts = 'XF:11BMB-ES{FC:1}F-I'
#FlowRate_SP = 'XF:11BMB-ES{FC:1}F:SP-SP'
#Mode_Sts = 'XF:11BMB-ES{FC:1}Mode:Opr-Sts'
#Mode_SP = 'XF:11BMB-ES{FC:1}Mode:Opr-Sel'
#ScaleFactor_SP = 'XF:11BMB-ES{FC:1}Val:ScaleFactor-SP'
#ScaleFactor_Sts = 'XF:11BMB-ES{FC:1}Val:ScaleFactor-RB'
#NominalRange_SP = 'XF:11BMB-ES{FC:1}F:FullRng-SP'
#NominalRange_Sts = 'XF:11BMB-ES{FC:1}F:FullRng-RB'

#FlowRate_Sts = 'XF:11BMB-ES{FC:1}F-I'
#FlowRate_SP = 'XF:11BMB-ES{FC:1}F:SP-SP'
#Mode_Sts = 'XF:11BMB-ES{FC:1}Mode:Opr-Sts'
#Mode_SP = 'XF:11BMB-ES{FC:1}Mode:Opr-Sel'
#ScaleFactor_SP = 'XF:11BMB-ES{FC:1}Val:ScaleFactor-SP'
#ScaleFactor_Sts = 'XF:11BMB-ES{FC:1}Val:ScaleFactor-RB'
#NominalRange_SP = 'XF:11BMB-ES{FC:1}F:FullRng-SP'
#NominalRange_Sts = 'XF:11BMB-ES{FC:1}F:FullRng-RB'
                  
class MassFlowControl(Device):
    
    def __init__(self):
        pass
        #self.FlowRate_Sts = 'XF:11BMB-ES{FC:1}F-I'
        #self.FlowRate_SP = 'XF:11BMB-ES{FC:1}F:SP-SP'
        #self.Mode_Sts = 'XF:11BMB-ES{FC:1}Mode:Opr-Sts'
        #self.Mode_SP = 'XF:11BMB-ES{FC:1}Mode:Opr-Sel'
        #self.ScaleFactor_SP = 'XF:11BMB-ES{FC:1}Val:ScaleFactor-SP'
        #self.ScaleFactor_Sts = 'XF:11BMB-ES{FC:1}Val:ScaleFactor-RB'
        #self.NominalRange_SP = 'XF:11BMB-ES{FC:1}F:FullRng-SP'
        #self.NominalRange_Sts = 'XF:11BMB-ES{FC:1}F:FullRng-RB'

    def setDevice(self, device='A1'):
        
        if device == 'A1':
            device_no = 1
        elif device == 'A2':
            device_no = 2 
        elif device == 'B1':
            device_no = 3
        elif device == 'B2':
            device_no = 4 
        else:
            print('The device is NOT valid.')
            device_no = 1
        print('Select Device {} in port {}'.format(device_no, device))
        self.FlowRate_Sts = 'XF:11BMB-ES{{FC:{}}}F-I'.format(device_no)
        self.FlowRate_SP = 'XF:11BMB-ES{{FC:{}}}F:SP-SP'.format(device_no)
        self.Mode_Sts = 'XF:11BMB-ES{{FC:{}}}Mode:Opr-Sts'.format(device_no)
        self.Mode_SP = 'XF:11BMB-ES{{FC:{}}}Mode:Opr-Sel'.format(device_no)
        self.ScaleFactor_SP = 'XF:11BMB-ES{{FC:{}}}Val:ScaleFactor-SP'.format(device_no)
        self.ScaleFactor_Sts = 'XF:11BMB-ES{{FC:{}}}Val:ScaleFactor-RB'.format(device_no)
        self.NominalRange_SP = 'XF:11BMB-ES{{FC:{}}}F:FullRng-SP'.format(device_no)
        self.NominalRange_Sts = 'XF:11BMB-ES{{FC:{}}}F:FullRng-RB'.format(device_no)
        
        return device_no
    
    def setFlow(self, device, rate, tolerence=1, verbosity=3):
        #set the setpoint of flow
        self.setDevice(device=device)
        caput(self.FlowRate_SP, rate)
        if verbosity >= 3:
            print('The flow rate has been set to {}'.format(rate))
            print('The current flow rate is {}'.format(caget(self.FlowRate_Sts)))
        return caget(self.FlowRate_Sts)
    
    def flow(self, device, verbosity=3):
        self.setDevice(device=device)
        if verbosity >= 3:
            print('The current flow rate is {}'.format(caget(self.FlowRate_Sts)))
        return caget(self.FlowRate_Sts)
        
    def setMode(self, device, mode, verbosity=3):
        #mode can be 0: OPEN, 1: Close, 2: SetPoint

        self.setDevice(device=device)

        if mode not in range(0, 3):
            return print('The input has to be 0: OPEN, 1: Close, 2: SetPoint')
        else:            
            caput(self.Mode_SP, mode)
            while self.mode(device=device) != mode:
                time.sleep(1)
                caput(self.Mode_SP, mode)
            if verbosity >=3:
                self.mode(device=device)
            return caget(self.Mode_Sts)
    
    def mode(self, device, verbosity=3):
        #OP mode: 0: OPEN, 1: Close, 2: SetPoint
        self.setDevice(device=device)
        if verbosity >= 3:
            self.readMode(device=device)
        return caget(self.Mode_Sts)
    
    def readMode(self, device, verbosity=3):
        #OP mode: 0: OPEN, 1: Close, 2: SetPoint
        self.setDevice(device=device)
        if caget(self.Mode_Sts)==0:
            return print('The current mode is {}'.format('ON'))
        if caget(self.Mode_Sts)==1:
            return print('The current mode is {}'.format('OFF'))
        if caget(self.Mode_Sts)==2:
            return print('The current mode is {}'.format('SETPOINT'))

    def scaleFactor(self, device, verbosity=3):
        #The scale factor depends on the gas. 
        #N2 and air are default as 1. Helium is 0.18.
        #More details are listed in the manual.
        self.setDevice(device=device)
        if verbosity >= 3:
            print('The scale factor is {}'.format(caget(self.ScaleFactor_Sts)))
        return caget(self.ScaleFactor_Sts)

    def setScaleFactor(self, device, val, verbosity=3):
        #Three modes: 0: OPEN, 1: Close, 2: SetPoint
        self.setDevice(device=device)

        caput(self.ScaleFactor_SP, val)
        while abs(self.scaleFactor() - val)>0.02:
            time.sleep(1)
            caput(self.ScaleFactor_SP, val)
        if verbosity >=3:
            self.scaleFactor()
        return caget(self.ScaleFactor_Sts)

    def deviceRange(self, device, verbosity=3):
        #Set the device by seting the Nominal Range. 
        #Model 201: 20SCCM. 
        self.setDevice(device=device)
        if verbosity >= 3:
            print('The range of this MFC device is up to {} SCCM.'.format(caget(self.NominalRange_Sts)))
        return caget(self.NominalRange_Sts)

    def setDeviceRange(self, device, val, verbosity=3):
        #Set the device by seting the Nominal Range. 
        #Model 201: 20SCCM. 
        self.setDevice(device=device)

        caput(self.NominalRange_SP, val)
        while abs(self.deviceRange() - val)>1:
            time.sleep(1)
            caput(self.NominalRange_SP, val)
        if verbosity >=3:
            self.deviceRange()
        return caget(self.NominalRange_Sts)            
            
        
    
ioL = ioLogik()
MFC = MassFlowControl()

