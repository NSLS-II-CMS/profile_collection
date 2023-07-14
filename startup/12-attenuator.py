"""
This code needs a bit of work before it it ready to go.
"""

from enum import Enum
from functools import reduce
from ophyd import Device, DeviceStatus, FormattedComponent, EpicsSignal, EpicsSignalRO, FormattedComponent
from ophyd.device import DynamicDeviceComponent


class StateEnum(Enum):
    In = True
    Out = False
    Unknown = None


class TernaryDevice(Device):
    """
    A general purpose ophyd device with set and reset signals, and a state signal
    with 3 posible signals.
    """

    set_cmd = FormattedComponent(EpicsSignal, "{self._set_name}")
    reset_cmd = FormattedComponent(EpicsSignal, "{self._reset_name}")
    state_rbv = FormattedComponent(EpicsSignalRO, "{self._state_name}", string=True)

    def __init__(
        self, *args, set_name, reset_name, state_name, state_enum, **kwargs
    ) -> None:
        self._state_enum = state_enum
        self._set_name = set_name
        self._reset_name = reset_name
        self._state_name = state_name
        self._state = None
        super().__init__(*args, **kwargs)

    def set(self, value=True):
        if value not in {True, False, 0, 1}:
            raise ValueError("value must be one of the following: True, False, 0, 1")

        target_value = bool(value)

        st = DeviceStatus(self)

        # If the device already has the requested state, return a finished status.
        if self._state == bool(value):
            st._finished()
            return st
        self._set_st = st

        def state_cb(value, timestamp, **kwargs):
            """
            Updates self._state and checks if the status should be marked as finished.
            """
            try:
                self._state = self._state_enum[value].value
            except KeyError:
                raise ValueError(f"self._state_enum does not contain value: {value}")
            if self._state == target_value:
                self._set_st = None
                self.state_rbv.clear_sub(state_cb)
                st._finished()

        # Subscribe the callback to the readback signal.
        # The callback will be called each time the PV value changes.
        self.state_rbv.subscribe(state_cb)

        # Write to the signal.
        if value:
            self.set_cmd.set(1)
        else:
            self.reset_cmd.set(1)
        return st

    def reset(self):
        self.set(False)

    def get(self):
        return self._state


def array_device_builder(components):
    class ArrayDeviceBase(Device):
        """
        An ophyd.Device that is an array of devices.

        The set method takes a list of values.
        the get method returns a list of values.
        Parameters
        ----------
        devices: iterable
            The array of ophyd devices.
        """
        ddc =  DynamicDeviceComponent(components)
        def __init__(self, *args, **kwargs):
            self.com
            super().__init__(*args,**kwargs)

        def set(self, values):
            if len(values) != len(self.devices):
                raise ValueError(
                    f"The number of values ({len(values)}) must match "
                    f"the number of devices ({len(self.devices)})"
                )

            # If the device already has the requested state, return a finished status.
            diff = [self.devices[i].get() != value for i, value in enumerate(values)]
            if not any(diff):
                return DeviceStatus(self)._finished()

            # Set the value of each device and return a union of the statuses.
            statuses = [self.devices[i].set(value) for i, value in enumerate(values)]
            st = reduce(lambda a, b: a & b, statuses)
            return st
    return ArrayDeviceBase


class CmsFilter(TernaryDevice):
    """
    This class is an example about how to create a TernaryDevice specialization
    for a specific implementation.
    """

    def __init__(self, index, *args, **kwargs):
        super().__init__(
            *args,
            set_name=f"XF:11BMB-OP{{Fltr:{index}}}Cmd:Opn-Cmd",
            reset_name=f"XF:11BMB-OP{{Fltr:{index}}}Cmd:Cls-Cmd",
            state_name=f"XF:11BMB-OP{{Fltr:{index}}}Pos-Sts",
            state_enum=StateEnum,
            **kwargs,
        )


#components = {f'c{i}': (CmsFilter, i, {}) for i in range(8)}
#attentuator = array_device_builder(components)(name='attenuator')
