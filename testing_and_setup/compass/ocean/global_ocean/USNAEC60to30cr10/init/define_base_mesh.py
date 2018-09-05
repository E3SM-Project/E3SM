#!/usr/bin/env python

import numpy as np
import coastal_tools as ct

def cellWidthVsLatLon():

  km = 1000.0

  params = ct.default_params
  cell_width,lon,lat = ct.coastal_refined_mesh(params)
  return cell_width/km,lon,lat
