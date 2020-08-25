#import time as ttime  # tea time
#from datetime import datetime
from ophyd import (ProsilicaDetector, SingleTrigger, TIFFPlugin,
                   ImagePlugin, StatsPlugin, DetectorBase, HDF5Plugin,
                   AreaDetector, EpicsSignal, EpicsSignalRO, ROIPlugin,
                   TransformPlugin, ProcessPlugin)
from ophyd.areadetector.cam import AreaDetectorCam
from ophyd.areadetector.base import ADComponent, EpicsSignalWithRBV
#from ophyd.areadetector.filestore_mixins import (FileStoreTIFFIterativeWrite,
#                                                 FileStoreHDF5IterativeWrite,
#                                                 FileStoreBase, new_short_uid)
from ophyd import Component as Cpt, Signal
from ophyd.utils import set_and_wait
#import filestore.api as fs


#class Elm(SingleTrigger, DetectorBase):
 #   pass




#class TIFFPluginWithFileStore(TIFFPlugin, FileStoreTIFFIterativeWrite):
#    pass


class StandardProsilica(SingleTrigger, ProsilicaDetector):
    # tiff = Cpt(TIFFPluginWithFileStore,
    #           suffix='TIFF1:',
    #           write_path_template='/XF11ID/data/')
    image = Cpt(ImagePlugin, 'image1:')
    stats1 = Cpt(StatsPlugin, 'Stats1:')
    stats2 = Cpt(StatsPlugin, 'Stats2:')
    stats3 = Cpt(StatsPlugin, 'Stats3:')
    stats4 = Cpt(StatsPlugin, 'Stats4:')
    stats5 = Cpt(StatsPlugin, 'Stats5:')
    trans1 = Cpt(TransformPlugin, 'Trans1:')
    roi1 = Cpt(ROIPlugin, 'ROI1:')
    roi2 = Cpt(ROIPlugin, 'ROI2:')
    roi3 = Cpt(ROIPlugin, 'ROI3:')
    roi4 = Cpt(ROIPlugin, 'ROI4:')
    proc1 = Cpt(ProcessPlugin, 'Proc1:')


#class StandardProsilicaWithTIFF(StandardProsilica):
#    tiff = Cpt(TIFFPluginWithFileStore,
#               suffix='TIFF1:',
#               write_path_template='/XF11ID/data/%Y/%m/%d/')



## This renaming should be reversed: no correspondance between CSS screens, PV names and ophyd....
#xray_eye1 = StandardProsilica('XF:11IDA-BI{Bpm:1-Cam:1}', name='xray_eye1')
#xray_eye2 = StandardProsilica('XF:11IDB-BI{Mon:1-Cam:1}', name='xray_eye2')
#xray_eye3 = StandardProsilica('XF:11IDB-BI{Cam:08}', name='xray_eye3')
#xray_eye1_writing = StandardProsilicaWithTIFF('XF:11IDA-BI{Bpm:1-Cam:1}', name='xray_eye1')
#xray_eye2_writing = StandardProsilicaWithTIFF('XF:11IDB-BI{Mon:1-Cam:1}', name='xray_eye2')
#xray_eye3_writing = StandardProsilicaWithTIFF('XF:11IDB-BI{Cam:08}', name='xray_eye3')
#fs1 = StandardProsilica('XF:11IDA-BI{FS:1-Cam:1}', name='fs1')
#fs2 = StandardProsilica('XF:11IDA-BI{FS:2-Cam:1}', name='fs2')
#fs_wbs = StandardProsilica('XF:11IDA-BI{BS:WB-Cam:1}', name='fs_wbs')
#dcm_cam = StandardProsilica('XF:11IDA-BI{Mono:DCM-Cam:1}', name='dcm_cam')
#fs_pbs = StandardProsilica('XF:11IDA-BI{BS:PB-Cam:1}', name='fs_pbs')
#elm = Elm('XF:11IDA-BI{AH401B}AH401B:',)

fs1 = StandardProsilica('XF:11BMA-BI{FS:1-Cam:1}', name='fs1')
fs2 = StandardProsilica('XF:11BMA-BI{FS:2-Cam:1}', name='fs2')
fs3 = StandardProsilica('XF:11BMB-BI{FS:3-Cam:1}', name='fs3')
fs4 = StandardProsilica('XF:11BMB-BI{FS:4-Cam:1}', name='fs4')
fs5 = StandardProsilica('XF:11BMB-BI{FS:Test-Cam:1}', name='fs5')

#all_standard_pros = [xray_eye1, xray_eye2, xray_eye3, xray_eye1_writing, xray_eye2_writing,
#                     xray_eye3_writing, fs1, fs2, fs_wbs, dcm_cam, fs_pbs]

all_standard_pros = [fs1, fs2, fs3, fs4, fs5]

for camera in all_standard_pros:
    camera.read_attrs = ['stats1', 'stats2','stats3','stats4','stats5']
    # camera.tiff.read_attrs = []  # leaving just the 'image'
    for stats_name in ['stats1', 'stats2','stats3','stats4','stats5']:
        stats_plugin = getattr(camera, stats_name)
        stats_plugin.read_attrs = ['total']
        camera.stage_sigs[stats_plugin.blocking_callbacks] = 1

    camera.stage_sigs[camera.roi1.blocking_callbacks] = 1
    camera.stage_sigs[camera.trans1.blocking_callbacks] = 1
    camera.stage_sigs[camera.cam.trigger_mode] = 'Fixed Rate'


#for camera in [xray_eye1_writing, xray_eye2_writing, xray_eye3_writing]:
#    camera.read_attrs.append('tiff')
#    camera.tiff.read_attrs = []


# Comment this out to suppress deluge of logging messages.
#import logging
#logging.basicConfig(level=logging.DEBUG)
#import ophyd.areadetector.filestore_mixins
#ophyd.areadetector.filestore_mixins.logger.setLevel(logging.DEBUG)
