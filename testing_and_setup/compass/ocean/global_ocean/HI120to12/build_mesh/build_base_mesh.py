#!/usr/bin/env python
"""
authors: Phillip J. Wolfram

This function specifies a high resolution patch for
Chris Jeffrey.
"""
import numpy as np
from mpas_tools.ocean import build_spherical_mesh


def cellWidthVsLatLon():
    """
    Create cell width array for this mesh on a regular latitude-longitude grid.
    Returns
    -------
       cellWidth : ndarray
            m x n array, entries are desired cell width in km
       lat : ndarray
            latitude, vector of length m, with entries between -90 and 90,
            degrees
       lon : ndarray
            longitude, vector of length n, with entries between -180 and 180,
            degrees
    """
    lat = np.arange(-90, 90.01, 1.0)
    lon = np.arange(-180, 180.01, 2.0)

    # in kms
    baseRes = 120.0
    highRes = 12.0
    latC = 20.0
    lonC = -155.0
    rad = 10.0

    theta = np.minimum(np.sqrt(((lat-latC)*(lat-latC))[:, np.newaxis] +
                               ((lon-lonC)*(lon-lonC))[np.newaxis, :])/rad, 1.0)

    cellWidth = (baseRes*theta + (1.0-theta)*highRes)*np.ones((lon.size,
                                                               lat.size))

    return cellWidth, lon, lat


def main():
    cellWidth, lon, lat = cellWidthVsLatLon()
    build_spherical_mesh(cellWidth, lon, lat, out_filename='base_mesh.nc',
                         do_inject_bathymetry=True)


if __name__ == '__main__':
    main()
