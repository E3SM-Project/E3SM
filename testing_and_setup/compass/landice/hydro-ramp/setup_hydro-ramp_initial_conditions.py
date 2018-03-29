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
   shiftx = unique_xs[len(unique_xs)/2]  # put origin in middle.
   xCell[:] = xCell[:] - shiftx
   xEdge[:] = xEdge[:] - shiftx
   xVertex[:] = xVertex[:] - shiftx



# thickness is a ramp
thickness[0,:] = (1000000.0-np.absolute(xCell[:])) * 0.001 + 0.0  # 0.0111 = 100 Pa/m
#thickness[0, thickness[0,:]<0.0 ]=0.0
thickness[0, np.absolute(xCell[:])>1000000.0]=0.0  # make a margin

# flat bed - make last cell floating in order to get N=0 lateral BC required
minH = np.unique(thickness[0,:])[1]  # 0.0 will the the smallset - get the next smallest
bedTopography[:] = -900.0/1000.0 * minH - 1.0  # subtract one extra meter to make sure we float here

# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

# melt
gridfile.variables['basalMeltInput'][:] = 2.0e-10 * 1000.0  # From Ian's email, 9/21/12
gridfile.variables['externalWaterInput'][:] = 2.0e-10 * 1000.0  # From Ian's email, 9/21/12

# velocity
gridfile.variables['uReconstructX'][:] = 1.0e-7  # doesn't matter because no sliding opening in the test case

# IC on thickness
gridfile.variables['waterThickness'][0,:] = (1000000.0 - np.absolute(xCell[:])) * 0.0 / 1000000.0
#gridfile.variables['waterThickness'][0,gridfile.variables['waterThickness'][0,:]<0.0] = 0.0

gridfile.close()

print 'Successfully added initial conditions to: ', options.filename

