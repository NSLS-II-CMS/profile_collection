print(f"Loading {__file__!r} ...")

##### Experimental shutters #####
# updated by RL, 20210901
# These shutters are controlled by sending a TTL pulse via Ecat controller.
from ophyd import Device
from ophyd.status import SubscriptionStatus

class TwoBladeShutter(Device):
    cmd = Cpt(EpicsSignal,"{Shutter}", kind="hinted", string=True)
    blade1_pos = Cpt(EpicsSignalRO, "{Psh_blade1}Pos", kind="hinted", string=True)
    blade2_pos = Cpt(EpicsSignalRO, "{Psh_blade2}Pos", kind="hinted", string=True)

    def set(self, value):
        if value not in {'OPEN', 'CLOSE'}:
            raise ValueError(f"Value: {value} must be either 'OPEN' or 'CLOSE")
        target = value

        def cb(*, value, **kwargs):
            return target == value
            
        st = SubscriptionStatus(self.blade1_pos, cb) & SubscriptionStatus(self.blade2_pos, cb)        
        self.cmd.set(value)
        return st
    
two_blade_shutter = TwoBladeShutter('XF:11BM-ES', name='two_blade_shutter')

#photonshutter_sts = EpicsSignal("XF:11BMA-PPS{PSh}Sts:FailOpn-Sts")
#photonshutter_open = EpicsSignal("XF:11BMA-PPS{PSh}Cmd:Opn-Cmd")
#photonshutter_cls = EpicsSignal("XF:11BMA-PPS{PSh}Cmd:Cls-Cmd")

def shutter_on(verbosity=3):
    yield from bps.mv(two_blade_shutter, 'OPEN')


def shutter_off(verbosity=3):
    yield from bps.mv(two_blade_shutter, 'CLOSE')


def shutter_state(verbosity=3):
    result = yield from bps.read(two_blade_shutter) 
    if result is None:
        return 1
    return int(all(sig['value'] == 'OPEN' for sig in result.values()))