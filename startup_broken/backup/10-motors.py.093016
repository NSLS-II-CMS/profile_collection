from ophyd import EpicsMotor, Device, Component as Cpt

slity = EpicsMotor('XF:11BMA-OP{Slt:0-Ax:T}Mtr', name='slity')


class Slits(Device):
    top = Cpt(EpicsMotor, '-Ax:T}Mtr')
    bottom = Cpt(EpicsMotor, '-Ax:B}Mtr')


slits = Slits('XF:11BMA-OP{Slt:0', name='slits') 
