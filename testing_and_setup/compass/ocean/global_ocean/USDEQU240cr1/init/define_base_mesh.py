#!/usr/bin/env python

import numpy as np
import mesh_definition_tools as mdt

def cellWidthVsLatLon():
  km = 1000.0
  
  params = mdt.params
  params["mesh_type"] = "QU"
  params["dx_max_global"] = 240.0*km
  params["region_box"] = mdt.Delaware_Bay
  params["dx_min_coastal"] = 1.0*km

  cell_width,lon,lat = mdt.coastal_refined_mesh(params)
  return cell_width/1000,lon,lat
