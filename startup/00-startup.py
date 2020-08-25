# import logging
# import caproto
# handler = logging.FileHandler('pilatus-trigger-log.txt')     
# from caproto._log import LogFormatter, color_log_format, log_date_format
# handler.setFormatter(                                                                                                                                                               
#     LogFormatter(color_log_format, datefmt=log_date_format))     
# caproto_log = logging.getLogger('caproto')                                                                                                                                                    
# caproto_log.handlers.clear()
# caproto_log.addHandler(handler)       
# logging.getLogger('caproto.ch').setLevel('DEBUG')
import nslsii
nslsii.configure_base(get_ipython().user_ns, 'cms')

from bluesky.magics import BlueskyMagics
from bluesky.preprocessors import pchain

# At the end of every run, verify that files were saved and
# print a confirmation message.
from bluesky.callbacks.broker import verify_files_saved
# RE.subscribe(post_run(verify_files_saved), 'stop')

from pyOlog.ophyd_tools import *

# Uncomment the following lines to turn on verbose messages for
# debugging.
# import logging
# ophyd.logger.setLevel(logging.DEBUG)
# logging.basicConfig(level=logging.DEBUG)

# Add a callback that prints scan IDs at the start of each scan.
# def print_scan_ids(name, start_doc):
#     print("Transient Scan ID: {0} @ {1}".format(start_doc['scan_id'],time.strftime("%Y/%m/%d %H:%M:%S")))
#     print("Persistent Unique Scan ID: '{0}'".format(start_doc['uid']))
# 
# RE.subscribe(print_scan_ids, 'start')

# - HACK #1 -  patch EpicsSignal.get to retry when timeouts happen stolen from HXN
import ophyd

def _epicssignal_get(self, *, as_string=None, connection_timeout=1.0, **kwargs):
    '''Get the readback value through an explicit call to EPICS
    Parameters
    ----------
    count : int, optional
        Explicitly limit count for array data
    as_string : bool, optional
        Get a string representation of the value, defaults to as_string
        from this signal, optional
    as_numpy : bool
        Use numpy array as the return type for array data.
    timeout : float, optional
        maximum time to wait for value to be received.
        (default = 0.5 + log10(count) seconds)
    use_monitor : bool, optional
        to use value from latest monitor callback or to make an
        explicit CA call for the value. (default: True)
    connection_timeout : float, optional
        If not already connected, allow up to `connection_timeout` seconds
        for the connection to complete.
    '''
    if as_string is None:
        as_string = self._string

    with self._lock:
        if not self._read_pv.connected:
            if not self._read_pv.wait_for_connection(connection_timeout):
                raise TimeoutError('Failed to connect to %s' %
                                   self._read_pv.pvname)

        ret = None
        attempts = 0
        max_attempts = 4
        while ret is None and attempts < max_attempts:
            attempts += 1
            ret = self._read_pv.get(as_string=as_string, **kwargs)
            if ret is None:
                print(f'*** PV GET TIMED OUT {self._read_pv.pvname} *** attempt #{attempts}/{max_attempts}')
        if ret is None:
            print(f'*** PV GET TIMED OUT {self._read_pv.pvname} *** return `None` as value :(')
            # TODO we really want to raise TimeoutError here, but that may cause more
            # issues in the codebase than we have the time to fix...
            # If this causes issues, remove it to keep the old functionality...
            raise TimeoutError('Failed to get %s after %d attempts' %
                               (self._read_pv.pvname, attempts))
        if attempts > 1:
            print(f'*** PV GET succeeded {self._read_pv.pvname} on attempt #{attempts}')

    if as_string:
        return ophyd.signal.waveform_to_string(ret)

    return ret


from ophyd import EpicsSignal
from ophyd import EpicsSignalRO
from ophyd.areadetector import EpicsSignalWithRBV

EpicsSignal.get = _epicssignal_get
EpicsSignalRO.get = _epicssignal_get
EpicsSignalWithRBV.get = _epicssignal_get