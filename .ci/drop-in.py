caput(beam.mono_bragg_pv, 1.03953)
pilatus2M.tiff.create_directory.set(-20)
sam = SampleGISAXS('test')
detselect(pilatus2M)

pilatus2M.cam.num_images.put(1)
