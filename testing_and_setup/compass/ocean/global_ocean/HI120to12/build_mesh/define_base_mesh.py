#!/usr/bin/env python
'''
name: define_base_mesh
authors: Phillip J. Wolfram

This function specifies a high resolution patch for 
Chris Jeffrey.

'''
import numpy as np

def cellWidthVsLatLon():
    lat = np.arange(-90, 90.01, 1.0)
    lon = np.arange(-180, 180.01, 2.0)

    km = 1000.0
    # in kms
    baseRes = 120.0
    highRes = 12.0
    latC = 20.0
    lonC = -155.0
    rad = 10.0

    theta = np.minimum(np.sqrt(((lat-latC)*(lat-latC))[:,np.newaxis] + ((lon-lonC)*(lon-lonC))[np.newaxis,:])/rad, 1.0)

    cellWidth = (baseRes*theta + (1.0-theta)*highRes)*np.ones((lon.size,lat.size))

    return cellWidth, lon, lat
