#!/usr/bin/env python
'''
This script plots results from MPAS-Ocean planar output.
'''
import numpy.ma as ma
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

varNames = ['temperature', 'salinity', 'tracer1', 'tracer2', 'tracer3']
nRow = len(varNames)
iTime = [0, 6, 12]
nCol = len(iTime)

fig = plt.gcf()
fig.set_size_inches(12.0, 10.0)

#ncfileMesh = Dataset('../base_mesh/planar_hex_mesh.nc','r')
# nx=ncfileMesh.getncattr('nx')
# ny=ncfileMesh.getncattr('ny')
# ncfileMesh.close()

ncfile = Dataset('output.nc', 'r')
ncfileIC = Dataset('../initial_state/initial_state.nc', 'r')
yMinkm = min(ncfile.variables['yCell']) / 1.0e3
yMaxkm = max(ncfile.variables['yCell']) / 1.0e3
zMin = 0.0
zMax = max(ncfile.variables['refBottomDepth'])

xtime = ncfile.variables['xtime']

# title at top:
# MPAS-O Test: Southern Ocean basin, 3000km x 4800m, cells:  40km x 100m
# Only Redi diffusion is on. All other tendencies are off. Nonlinear EOS.
# slope: 0.01

titleTxt = [', initial', ', 6 months', ', 1 year']
for iRow in range(nRow):
    var = ncfile.variables[varNames[iRow]]
    for iCol in range(nCol):
        plt.subplot(nRow, nCol, iRow * nCol + iCol + 1)
        varSliced = np.transpose(var[iTime[iCol], ::4, :])
        varMasked = ma.masked_where(varSliced < -1.0e33, varSliced)
        ax = plt.imshow(varMasked)  # ,extent=[yMinkm,yMaxkm,zMin,zMax])
        plt.axis('off')
        plt.title(varNames[iRow] + titleTxt[iCol])
        plt.jet()
        if iRow == nRow - 1:
            plt.xlabel('x, km')
        if iCol == 0:
            plt.ylabel('z, m')
        plt.colorbar()

ncfile.close()
ncfileIC.close()
plt.savefig('output.png')
