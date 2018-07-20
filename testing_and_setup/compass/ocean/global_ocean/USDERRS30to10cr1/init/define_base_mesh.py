#!/usr/bin/env python

import numpy as np
import mesh_definition_tools as mdt

def cellWidthVsLatLon():

  params = mdt.params
  params["mesh_type"] = "RRS"
  params["region_box"] = mdt.Delaware_Bay
  params["dx_min_coastal"] = 1000

  cell_width,lon,lat = mdt.coastal_refined_mesh(params)
  return cell_width/1000,lon,lat
