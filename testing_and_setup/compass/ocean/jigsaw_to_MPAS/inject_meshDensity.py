#!/usr/bin/env python
# Simple script to inject mesh density onto a mesh
# example usage:
#   ./inject_meshDensity.py cellWidthVsLatLon.mat base_mesh.nc
# where:
#   cellWidthVsLatLon.mat is a matlab workspace file with cell width
#   base_mesh.nc is the mpas netcdf file where meshDensity is added
# Mark Petersen, 7/24/2018

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import netCDF4 as nc4
import scipy.io as sio
dtor = np.pi / 180.0
rtod = 180.0 / np.pi

if __name__ == "__main__":
    import sys

    print 'Read cell width field from Matlab regular grid...'
    matData = sio.loadmat(sys.argv[1])
    # Add extra column in longitude to interpolate over the Date Line
    cellWidth = np.concatenate(
        (matData['cellWidth'], matData['cellWidth'][:, 0:1]), axis=1)
    LonPos = np.concatenate(
        (matData['lon'].T, matData['lon'].T[0:1] + 360)) * dtor
    LatPos = matData['lat'].T * dtor
    # set max lat position to be exactly at North Pole to avoid interpolation
    # errors
    LatPos[np.argmax(LatPos)] = np.pi / 2.0
    minCellWidth = cellWidth.min()
    meshDensityVsLatLon = (minCellWidth / cellWidth)**4
    print '  minimum cell width in grid definition:', minCellWidth
    print '  maximum cell width in grid definition:', cellWidth.max()

    LON, LAT = np.meshgrid(LonPos, LatPos)

    print 'Open unstructured MPAS mesh file...'
    ds = nc4.Dataset(sys.argv[2], 'r+')
    meshDensity = ds.variables['meshDensity'][:]

    print 'Preparing interpolation of meshDensity from lat/lon to mesh...'
    meshDensityInterp = interpolate.LinearNDInterpolator(
        np.vstack((LAT.ravel(), LON.ravel())).T, meshDensityVsLatLon.ravel())

    print 'Interpolating and writing meshDensity...'
    ds.variables['meshDensity'][:] = meshDensityInterp(
        np.vstack(
            (ds.variables['latCell'][:],
             np.mod(
                ds.variables['lonCell'][:] +
                np.pi,
                2 *
                np.pi) -
                np.pi)).T)

    ds.close()
