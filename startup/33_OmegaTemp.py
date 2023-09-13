
print(f"Loading {__file__!r} ...")

from epics import caput, caget


class OmegaTempCntl78000(Device):
    temperature_current = Cpt(EpicsSignal, '{LA:Omega}_T_RBV')
    temperature_setpoint = Cpt(EpicsSignal, '{LA:Omega}_SP')

    def setTemperature(self, temperature):
        return self.temperature_setpoint.put(temperature)

tempCntlOmega = OmegaTempCntl78000("XF:11BM", name='omega')

