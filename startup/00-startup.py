print(f"Loading {__file__!r} ...")

import nslsii
import os
from bluesky.magics import BlueskyMagics
from bluesky.preprocessors import pchain
from bluesky.utils import PersistentDict
from pyOlog.ophyd_tools import *

# At the end of every run, verify that files were saved and
# print a confirmation message.
from bluesky.callbacks.broker import verify_files_saved


# Added this variable temporarily to bypass some code that doesn't run without the beamline.
# This can be set when starting bsui like this: `BS_MODE=test bsui`

testing = os.environ.get("BS_TESTING", False) == "True"
PTA_MODE = os.environ.get("BS_PTA", False) == "True"


# This is the definition for cms object.
CMS_BEAMLINE_MODE = "GISAXS"
if CMS_BEAMLINE_MODE not in ["Standard", "XR", "GISAXS"]:
    raise ValueError("CMS_BEAMLINE MODE must be one of the following values: Standard, XR, GISAXS")
print(f"{CMS_BEAMLINE_MODE =}")

# This is the definition for the sample stage.
CMS_BEAMLINE_STAGE = "default"
if CMS_BEAMLINE_STAGE not in ["default", "open_MAXS", "BigHuber", "testing"]:
    raise ValueError(
        "CMS_BEAMLINE MODE must be one of the following values: default, open_MAXS, BigHuber, testing"
    )
print(f"{CMS_BEAMLINE_STAGE =}")


if testing:
    nslsii.configure_base(get_ipython().user_ns, "temp", publish_documents_with_kafka=False)
else:
    nslsii.configure_base(get_ipython().user_ns, "cms", publish_documents_with_kafka=True)

# RE.subscribe(post_run(verify_files_saved), 'stop')


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

# runengine_metadata_dir = appdirs.user_data_dir(appname="bluesky") / Path("runengine-metadata")
if testing:
    runengine_metadata_dir = "/tmp/runingine-metadata"
else:
    runengine_metadata_dir = "/nsls2/data/cms/shared/config/runengine-metadata"

# PersistentDict will create the directory if it does not exist
RE.md = PersistentDict(runengine_metadata_dir)

# The following plan stubs are automatically imported in global namespace by 'nslsii.configure_base',
# but have signatures that are not compatible with the Queue Server. They should not exist in the global
# namespace, but can be accessed as 'bps.one_1d_step' etc. from other plans.
del one_1d_step, one_nd_step, one_shot
