#!/usr/bin/env python

import numpy as np
import mesh_definition_tools as mdt

def cellWidthVsLatLon():
  params = mdt.params
  cell_width,lon,lat = mdt.coastal_refined_mesh(params)
  return cell_width/1000,lon,lat
