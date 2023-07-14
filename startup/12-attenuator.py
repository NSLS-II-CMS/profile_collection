"""
This code is ready to be tested.

Here we create an ArrayDevice for the attentuator.

This allows us to do:

attenuator.set([1,0,1])

and:

attenuator.get() will return an array with all of the filter states.
"""

from enum import Enum
from functools import reduce
from ophyd import Device, DeviceStatus, FormattedComponent, EpicsSignal, EpicsSignalRO, FormattedComponent
from ophyd.device import DynamicDeviceComponent


class TernaryDevice(Device):
    """
    A general purpose ophyd device with set and reset signals, and a state signal
    with 3 posible signals.
    """

    set_cmd = FormattedComponent(EpicsSignal, "{self._set_name}")
    reset_cmd = FormattedComponent(EpicsSignal, "{self._reset_name}")
    state_rbv = FormattedComponent(EpicsSignalRO, "{self._state_name}", kind='hinted')

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
        return self.state_rbv.get()


def ArrayDevice(components, *args, **kwargs):
    """
    A function, that behaves like a class init, that dynamically creates an
    ArrayDevice class. This is needed to set class attributes before the init.
    Adding devices in the init can subvert important ophyd code that
    manages sub devices.
    """

    class _ArrayDeviceBase(Device):
        """
        An ophyd.Device that is an array of devices.

        The set method takes a list of values.
        the get method returns a list of values.
        Parameters
        ----------
        devices: iterable
            The array of ophyd devices.
        """
        def set(self, values):
            if len(values) != len(self.component_names):
                raise ValueError(
                    f"The number of values ({len(values)}) must match "
                    f"the number of devices ({len(self.devices)})"
                )

            # If the device already has the requested state, return a finished status.
            diff = [getattr(self, self.devices[i]).get() != value for i, value in enumerate(values)]
            if not any(diff):
                return DeviceStatus(self)._finished()

            # Set the value of each device and return a union of the statuses.
            statuses = [getattr(self, self.devices[i]).set(value) for i, value in enumerate(values)]
            st = reduce(lambda a, b: a & b, statuses)
            return st

        def get(self):
            return [getattr(self, device).get() for device in self.devices]

    #types = {component.cls for component in components}
    #if len(types) != 1:
    #    raise TypeError("All components must have the same type")

    clsdict = OrderedDict(
        {
            **{'__doc__': "ArrayDevice",
               '_default_read_attrs': None,
               '_default_configuration_attrs': None},
            **components,
            'devices': list(components.keys())
        }
    )

    _ArrayDevice = type('ArrayDevice', (_ArrayDeviceBase,), clsdict)
    return _ArrayDevice(*args, **kwargs)


class CmsEnum(Enum):
    In = True
    Out = False
    Unknown = None


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
            state_enum=CmsEnum,
            **kwargs,
        )

# TODO: Uncomment these lines when ready to test.
#filters = {f'filter{i}': FormattedComponent(CmsFilter, f"{i}") for i in range(8)}
#attenuator = ArrayDevice(filters, name="attenuator")
