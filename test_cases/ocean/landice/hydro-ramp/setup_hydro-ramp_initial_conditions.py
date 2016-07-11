#!/usr/bin/env python
# Generate initial conditions for hydro-margin land ice test case

import sys
from netCDF4 import Dataset as NetCDFFile
import numpy as np

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file to setup dome", metavar="FILE")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'


# Open the file, get needed dimensions
gridfile = NetCDFFile(options.filename,'r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
# Get variables
xCell = gridfile.variables['xCell']
yCell = gridfile.variables['yCell']
xEdge = gridfile.variables['xEdge']
yEdge = gridfile.variables['yEdge']
xVertex = gridfile.variables['xVertex']
yVertex = gridfile.variables['yVertex']
thickness = gridfile.variables['thickness']
bedTopography = gridfile.variables['bedTopography']
layerThicknessFractions = gridfile.variables['layerThicknessFractions']
SMB = gridfile.variables['sfcMassBal']


unique_xs=np.array(sorted(list(set(xCell[:]))))
if min(unique_xs) >= 0.0:
   shiftx = unique_xs[-5]  # put margin at 0, a few cells from the right edge
   xCell[:] = xCell[:] - shiftx
   xEdge[:] = xEdge[:] - shiftx
   xVertex[:] = xVertex[:] - shiftx


ind=xCell[:]<-1000000.0  # locations west of area of interest

# thickness is a ramp
thickness[0,:] = -1.0 * xCell[:] * 0.001 + 0.0  # 0.0111 = 100 Pa/m
thickness[0,thickness[0,:]<0.0]=0.0
thickness[0,ind] = 1000000.0 * 0.001  # make flat
thickness[0,xCell[:]<-1040000.0]=0.0  # make a margin

# flat bed
bedTopography[:] = 0.0

# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

# melt
#gridfile.variables['meltInput'][:] = 2.0e-9  + 0.06/3.0e5 + (1.0e-7 * 10.0**5)   # From Ian Hewitt email
gridfile.variables['meltInput'][:] = 2.0e-7
gridfile.variables['meltInput'][0,ind] = 0.0  # remove input to left of area of interest

# velocity
gridfile.variables['uReconstructX'][:] = 1.0e-7
gridfile.variables['uReconstructX'][0,ind,:] = 0.0
gridfile.variables['uReconstructX'][0,thickness[0,:]==0.0,:] = 0.0

# IC on thickness
gridfile.variables['waterThickness'][0,:] = (1000000.0 + xCell[:]) * 0.4 / 1000000.0
gridfile.variables['waterThickness'][0,gridfile.variables['waterThickness'][0,:]<0.0] = 0.0

gridfile.close()

print 'Successfully added initial conditions to: ', options.filename

