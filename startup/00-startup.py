import nslsii
from bluesky.magics import BlueskyMagics
from bluesky.preprocessors import pchain
from pyOlog.ophyd_tools import *

# At the end of every run, verify that files were saved and
# print a confirmation message.
from bluesky.callbacks.broker import verify_files_saved


nslsii.configure_base(get_ipython().user_ns, "cms", publish_documents_with_kafka=True)


# RE.subscribe(post_run(verify_files_saved), 'stop')


# Added this variable temporarily to bypass some code that doesn't run without the beamline.
testing = True

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
try:
    from bluesky.utils import PersistentDict
except ImportError:
    import msgpack
    import msgpack_numpy
    import zict

    class PersistentDict(zict.Func):
        """
        A MutableMapping which syncs it contents to disk.
        The contents are stored as msgpack-serialized files, with one file per item
        in the mapping.
        Note that when an item is *mutated* it is not immediately synced:
        >>> d['sample'] = {"color": "red"}  # immediately synced
        >>> d['sample']['shape'] = 'bar'  # not immediately synced
        but that the full contents are synced to disk when the PersistentDict
        instance is garbage collected.
        """

        def __init__(self, directory):
            self._directory = directory
            self._file = zict.File(directory)
            self._cache = {}
            super().__init__(self._dump, self._load, self._file)
            self.reload()

            # Similar to flush() or _do_update(), but without reference to self
            # to avoid circular reference preventing collection.
            # NOTE: This still doesn't guarantee call on delete or gc.collect()!
            #       Explicitly call flush() if immediate write to disk required.
            def finalize(zfile, cache, dump):
                zfile.update((k, dump(v)) for k, v in cache.items())

            import weakref

            self._finalizer = weakref.finalize(self, finalize, self._file, self._cache, PersistentDict._dump)

        @property
        def directory(self):
            return self._directory

        def __setitem__(self, key, value):
            self._cache[key] = value
            super().__setitem__(key, value)

        def __getitem__(self, key):
            return self._cache[key]

        def __delitem__(self, key):
            del self._cache[key]
            super().__delitem__(key)

        def __repr__(self):
            return f"<{self.__class__.__name__} {dict(self)!r}>"

        @staticmethod
        def _dump(obj):
            "Encode as msgpack using numpy-aware encoder."
            # See https://github.com/msgpack/msgpack-python#string-and-binary-type
            # for more on use_bin_type.
            return msgpack.packb(obj, default=msgpack_numpy.encode, use_bin_type=True)

        @staticmethod
        def _load(file):
            return msgpack.unpackb(file, object_hook=msgpack_numpy.decode, raw=False)

        def flush(self):
            """Force a write of the current state to disk"""
            for k, v in self.items():
                super().__setitem__(k, v)

        def reload(self):
            """Force a reload from disk, overwriting current cache"""
            self._cache = dict(super().items())


# runengine_metadata_dir = appdirs.user_data_dir(appname="bluesky") / Path("runengine-metadata")
runengine_metadata_dir = "/nsls2/data/cms/shared/config/runengine-metadata"


# PersistentDict will create the directory if it does not exist
RE.md = PersistentDict(runengine_metadata_dir)

print("a new version of bsui")
print("sth is happening")

# The following plan stubs are automatically imported in global namespace by 'nslsii.configure_base',
# but have signatures that are not compatible with the Queue Server. They should not exist in the global
# namespace, but can be accessed as 'bps.one_1d_step' etc. from other plans.
del one_1d_step, one_nd_step, one_shot
