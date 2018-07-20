#/usr/bin/env python

import mesh_definition_tools as mdt

def cellWidthVsLatLon():
    params = mdt.params
    lon,lat,cell_width = mdt.create_background_mesh(params["grd_box"],params["ddeg"],params["mesh_type"],params["dx_max"])
    cell_width = cell_width/1000
  
    return cell_width,lon,lat
