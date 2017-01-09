#!/usr/bin/env python
# Simple script to inject bathymetry onto a mesh
# Phillip Wolfram, 01/19/2018

import matplotlib.pyplot as plt
from open_msh import readmsh
import numpy as np
from scipy import interpolate
import netCDF4 as nc4
dtor = np.pi/180.0
rtod = 180.0/np.pi

if __name__ == "__main__":
    import sys

    topo = readmsh('topo.msh')
    xpos = topo['COORD1']*dtor
    ypos = topo['COORD2']*dtor
    zlev = np.reshape(topo['VALUE'], (len(ypos), len(xpos)))

    Y, X = np.meshgrid(ypos, xpos)

    ds = nc4.Dataset(sys.argv[1],'r+')
    ds.createVariable('bathymetry','f8',('nCells'))
    ds.createVariable('cullCell','i',('nCells'))
    bathy = interpolate.LinearNDInterpolator(np.vstack((X.ravel(), Y.ravel())).T, zlev.ravel())
    ds.variables['bathymetry'][:] = bathy(np.vstack((np.mod(ds.variables['lonCell'][:] + np.pi, 2*np.pi)-np.pi,
                                          ds.variables['latCell'][:])).T)
    ds.variables['cullCell'][:] = ds.variables['bathymetry'][:] > 20.0

    ds.close()


