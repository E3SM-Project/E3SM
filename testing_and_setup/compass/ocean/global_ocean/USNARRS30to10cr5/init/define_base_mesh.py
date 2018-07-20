#!/usr/bin/env python

import numpy as np
import mesh_definition_tools as mdt

def cellWidthVsLatLon():
  params = mdt.params
  params["mesh_type"] = "RRS"
  params["dx_min_coastal"] = 5000.0
  cell_width,lon,lat = mdt.coastal_refined_mesh(params)
  return cell_width/1000.0,lon,lat
