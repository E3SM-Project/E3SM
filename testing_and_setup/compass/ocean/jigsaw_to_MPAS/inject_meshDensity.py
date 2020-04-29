#!/usr/bin/env python
# Simple script to inject mesh density onto a mesh
# example usage:
#   ./inject_meshDensity.py cellWidthVsLatLon.nc base_mesh.nc
# where:
#   cellWidthVsLatLon.nc is a netcdf file with cellWidth
#   base_mesh.nc is the mpas netcdf file where meshDensity is added
# Mark Petersen, 7/24/2018

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
from scipy import interpolate
import netCDF4 as nc4
import sys


def inject_meshDensity(cw_filename, mesh_filename, on_sphere=True):
    print('Read cell width field from nc file regular grid...')
    ds = nc4.Dataset(cw_filename,'r')
    cellWidth = ds.variables['cellWidth'][:]
    if on_sphere:
       lon = ds.variables['lon'][:]
       lat = ds.variables['lat'][:]
    else:
       x = ds.variables['x'][:]
       y = ds.variables['y'][:]
    ds.close()

    if on_sphere:
       # Add extra column in longitude to interpolate over the Date Line
       cellWidth = np.concatenate(
           (cellWidth, cellWidth[:, 0:1]), axis=1)
       LonPos = np.deg2rad(np.concatenate(
           (lon.T, lon.T[0:1] + 360)))
       LatPos = np.deg2rad(lat.T)
       # set max lat position to be exactly at North Pole to avoid interpolation
       # errors
       LatPos[np.argmax(LatPos)] = np.pi / 2.0
    minCellWidth = cellWidth.min()
    meshDensityVsXY = (minCellWidth / cellWidth)**4
    print('  minimum cell width in grid definition: {0:.0f} km'.format(minCellWidth/1000.))
    print('  maximum cell width in grid definition: {0:.0f} km'.format(cellWidth.max()/1000.))

    if on_sphere:
       X, Y = np.meshgrid(LonPos, LatPos)
    else:
       X, Y = np.meshgrid(x, y)

    print('Open unstructured MPAS mesh file...')
    ds = nc4.Dataset(mesh_filename, 'r+')
    meshDensity = ds.variables['meshDensity']

    print('Preparing interpolation of meshDensity from native coordinates to mesh...')
    meshDensityInterp = interpolate.LinearNDInterpolator(
        np.vstack((X.ravel(), Y.ravel())).T, meshDensityVsXY.ravel())

    print('Interpolating and writing meshDensity...')
    if on_sphere:
       meshDensity[:] = meshDensityInterp(
           np.vstack(
            (np.mod(
                ds.variables['lonCell'][:] +
                np.pi,
                2 *
                np.pi) -
                np.pi,
             ds.variables['latCell'][:])).T)
    else:
       meshDensity[:] = meshDensityInterp(
           np.vstack(
            (ds.variables['xCell'][:],
             ds.variables['yCell'][:])
            ).T)

    ds.close()


if __name__ == "__main__":

    inject_meshDensity(cw_filename=sys.argv[1], mesh_filename=sys.argv[2])
