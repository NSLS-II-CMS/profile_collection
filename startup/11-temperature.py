print(f"Loading {__file__!r} ...")

from ophyd import PVPositionerPC
from ophyd import FormattedComponent
from ophyd import EpicsSignal, EpicsSignalRO
from ophyd import Component as Cpt


class Lakeshore(PVPositionerPC):
    """
    Ophyd device for Lakeshore that uses put completion.
    PVPositionerPC does not require a done signal like PVPositioner,
    instead it uses the setpoint put_completion.

    Example
    -------
    lakeshore = Lakeshore('Lakeshore', output="1", chan="A", settle_time=5)
    lakeshore.set(100).wait()
    This will wait for the ramp to be completed and also wait for
    the settle_time.

    lakeshore.get() to see all of the pv values.
    """

    feedback = FormattedComponent(EpicsSignalRO, "{self.prefix}-Chan:{self._chan}}}T:C-I")
    output = FormattedComponent(EpicsSignalRO, "{self.prefix}-Out:{self._output}}}Out-I")
    setpoint = FormattedComponent(EpicsSignal, "{self.prefix}-Out:{self._output}}}T-SP", put_complete=True)
    ramp_rate = FormattedComponent(EpicsSignal, "{self.prefix}-Out:{self._output}}}Val:Ramp-SP")

    def __init__(self, *args, chan, output, **kwargs) -> None:
        self._chan = chan
        self._output = output
        super().__init__(*args, **kwargs)


lakeshore = Lakeshore("XF:11BM-ES{Env:01", output="1", chan="A", settle_time=5)
