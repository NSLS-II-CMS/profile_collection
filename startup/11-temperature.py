from ophyd import PVPositionerPC


class Lakeshore(PVPositionerPC):
    """
    Ophyd device for Lakeshore that uses put completion.
    PVPositionerPC does not require a done signal like PVPositioner,
    instead it uses the setpoint put_completion.
    
    Example
    -------
    ls = Lakeshore('Lakeshore', name='Lakeshore', settle_time=5)
    ls.set(100).wait()
    This will wait for the ramp to be completed and also wait for
    the settle_time.
    """

    feedback = Cpt(EpicsSignalRO, ":feedback")
    output = Cpt(EpicsSignalRO, ":output")
    setpoint = Cpt(EpicsSignal, ":setpoint", put_complete=True)
    ramp_rate = Cpt(EpicsSignal, ":ramp_rate")
