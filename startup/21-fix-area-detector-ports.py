for fs in [fs1, fs2, fs3, fs4, fs5]:
    G, port_dict = fs.get_asyn_digraph()
    cam_port = fs.cam.port_name.get()
    for v in port_dict.values():
        try:
            if v.nd_array_port.get() == "CAM":
                v.nd_array_port.put(cam_port)
        except AttributeError:
            pass
    fs.validate_asyn_ports()
