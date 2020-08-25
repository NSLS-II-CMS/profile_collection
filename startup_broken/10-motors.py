from ophyd import EpicsMotor, Device, Component as Cpt

#slity = EpicsMotor('XF:11BMA-OP{Slt:0-Ax:T}Mtr', name='slity')


#class Slits(Device):
#    top = Cpt(EpicsMotor, '-Ax:T}Mtr')
#    bottom = Cpt(EpicsMotor, '-Ax:B}Mtr')


#slits = Slits('XF:11BMA-OP{Slt:0', name='slits') 

###################################################################################
#above as found on 19 Oct 2016, then commented out
#below added by MF in Oct 2016
###################################################################################


########## motor classes ##########
class MotorCenterAndGap(Device):
    "Center and gap using Epics Motor records"
    xc = Cpt(EpicsMotor, '-Ax:XC}Mtr')
    yc = Cpt(EpicsMotor, '-Ax:YC}Mtr')
    xg = Cpt(EpicsMotor, '-Ax:XG}Mtr')
    yg = Cpt(EpicsMotor, '-Ax:YG}Mtr')

class Blades(Device):
    "Actual T/B/O/I and virtual center/gap using Epics Motor records"
    tp = Cpt(EpicsMotor, '-Ax:T}Mtr')
    bt = Cpt(EpicsMotor, '-Ax:B}Mtr')
    ob = Cpt(EpicsMotor, '-Ax:O}Mtr')
    ib = Cpt(EpicsMotor, '-Ax:I}Mtr')
    xc = Cpt(EpicsMotor, '-Ax:XCtr}Mtr')
    yc = Cpt(EpicsMotor, '-Ax:YCtr}Mtr')
    xg = Cpt(EpicsMotor, '-Ax:XGap}Mtr')
    yg = Cpt(EpicsMotor, '-Ax:YGap}Mtr')

#class MotorSlits(Blades, MotorCenterAndGap):
#    "combine t b i o and xc yc xg yg"
#    pass

#class VirtualMotorSlits(Blades, VirtualMotorCenterAndGap):
#    "combine t b i o and xc yc xg yg"
#    pass


########## FOE motors ##########
## stages for monochromator (FOE)
mono_bragg = EpicsMotor('XF:11BMA-OP{Mono:DMM-Ax:Bragg}Mtr', name='mono_bragg')
mono_pitch2 = EpicsMotor('XF:11BMA-OP{Mono:DMM-Ax:P2}Mtr', name='mono_pitch2')
mono_roll2 = EpicsMotor('XF:11BMA-OP{Mono:DMM-Ax:R2}Mtr', name='mono_roll2')
mono_perp2 = EpicsMotor('XF:11BMA-OP{Mono:DMM-Ax:Y2}Mtr', name='mono_perp2')

## stages for toroidal mirror (FOE)
mir_usx = EpicsMotor('XF:11BMA-OP{Mir:Tor-Ax:XU}Mtr', name='mir_usx')
mir_dsx = EpicsMotor('XF:11BMA-OP{Mir:Tor-Ax:XD}Mtr', name='mir_dsx')
mir_usy = EpicsMotor('XF:11BMA-OP{Mir:Tor-Ax:YU}Mtr', name='mir_usy')
mir_dsyi = EpicsMotor('XF:11BMA-OP{Mir:Tor-Ax:YDI}Mtr', name='mir_dsyi')
mir_dsyo = EpicsMotor('XF:11BMA-OP{Mir:Tor-Ax:YDO}Mtr', name='mir_dsyo')
mir_bend = EpicsMotor('XF:11BMA-OP{Mir:Tor-Ax:UB}Mtr', name='mir_bend')


########## FOE slits ##########
## FMB Oxford slits -- usage: s0.tp, s0.bt, s0.ob, s0.ib, s0.xc, s0.xg, s0.yc, s0.yg
s0 = Blades('XF:11BMA-OP{Slt:0', name='s0') 	


########## Endstation slits #########		
## jj slits -- usage: s*.xc, s*.xg, s*.yc, s*.yg
s1 = MotorCenterAndGap('XF:11BMB-OP{Slt:1', name='s1')
s2 = MotorCenterAndGap('XF:11BMB-OP{Slt:2', name='s2') 
s3 = MotorCenterAndGap('XF:11BMB-OP{Slt:3', name='s3') 
s4 = MotorCenterAndGap('XF:11BMB-OP{Slt:4', name='s4') 
s5 = MotorCenterAndGap('XF:11BMB-OP{Slt:5', name='s5') 


########## Endstation motors ##########
## stages for Endstation diagnostics
bim3y = EpicsMotor('XF:11BMB-BI{IM:3-Ax:Y}Mtr', name='bim3y')
fs3y = EpicsMotor('XF:11BMB-BI{FS:3-Ax:Y}Mtr', name='fs3y')
bim4y = EpicsMotor('XF:11BMB-BI{IM:4-Ax:Y}Mtr', name='bim4y')
bim5y = EpicsMotor('XF:11BMB-BI{IM:5-Ax:Y}Mtr', name='bim5y')

## stages for sample positioning
smx = EpicsMotor('XF:11BMB-ES{Chm:Smpl-Ax:X}Mtr', name='smx')
smy = EpicsMotor('XF:11BMB-ES{Chm:Smpl-Ax:Z}Mtr', name='smy')
sth = EpicsMotor('XF:11BMB-ES{Chm:Smpl-Ax:theta}Mtr', name='sth')
schi = EpicsMotor('XF:11BMB-ES{Chm:Smpl-Ax:chi}Mtr', name='schi')
sphi = EpicsMotor('XF:11BMB-ES{Chm:Smpl-Ax:phi}Mtr', name='sphi')
srot = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Srot}Mtr', name='srot')
strans = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Strans}Mtr', name='strans')
strans2 = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Strans2}Mtr', name='strans2')
stilt = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Stilt}Mtr', name='stilt')
stilt2 = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Stilt2}Mtr', name='stilt2')


## stages for on-axis sample camera mirror/lens
camx = EpicsMotor('XF:11BMB-ES{Cam:OnAxis-Ax:X1}Mtr', name='camx')
camy = EpicsMotor('XF:11BMB-ES{Cam:OnAxis-Ax:Y1}Mtr', name='camy')

## stages for off-axis sample camera
cam2x = EpicsMotor('XF:11BMB-ES{Cam:OnAxis-Ax:X2}Mtr', name='cam2x')
cam2z = EpicsMotor('XF:11BMB-ES{Cam:OnAxis-Ax:Y2}Mtr', name='cam2z')

## stages for sample exchanger
armz = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Z}Mtr', name='armz')
armx = EpicsMotor('XF:11BMB-ES{SM:1-Ax:X}Mtr', name='armx')
armphi = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Yaw}Mtr', name='armphi')
army = EpicsMotor('XF:11BMB-ES{SM:1-Ax:Y}Mtr', name='army')
armr = EpicsMotor('XF:11BMB-ES{SM:1-Ax:ArmR}Mtr', name='armr')

## stages for detectors
DETx = EpicsMotor('XF:11BMB-ES{Det:Stg-Ax:X}Mtr', name='DETx')
DETy =  EpicsMotor('XF:11BMB-ES{Det:Stg-Ax:Y}Mtr', name='DETy')
WAXSx = EpicsMotor('XF:11BMB-ES{Det:WAXS-Ax:X}Mtr', name='WAXSx')

SAXSx = EpicsMotor('XF:11BMB-ES{Det:SAXS-Ax:X}Mtr', name='SAXSx')
SAXSy = EpicsMotor('XF:11BMB-ES{Det:SAXS-Ax:Y}Mtr', name='SAXSy')

## stages for beamstops
bsx = EpicsMotor('XF:11BMB-ES{BS:SAXS-Ax:X}Mtr', name='bsx')
bsy = EpicsMotor('XF:11BMB-ES{BS:SAXS-Ax:Y}Mtr', name='bsy')
bsphi = EpicsMotor('XF:11BMB-ES{BS:SAXS-Ax:phi}Mtr', name='bsphi')


