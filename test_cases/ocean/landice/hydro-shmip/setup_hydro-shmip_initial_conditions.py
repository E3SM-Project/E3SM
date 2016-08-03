#!/usr/bin/env python
# Generate initial conditions for hydro-SHMIP land ice test cases
# Details here: http://shmip.bitbucket.org/

from netCDF4 import Dataset as NetCDFFile
import numpy as np

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file to setup", metavar="FILE")
parser.add_option("-t", "--test", dest="test", type='string', help="test type to set up, A-F", metavar="LETTER")
parser.add_option("-n", "--number", dest="number", type='int', help="test variant to set up, 1-X, depending on test type", metavar="NUMBER")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'

# Setup dictionaries of parameter values for each experiment
a_params = {1:7.93e-11,  2:1.59e-09, 3:5.79e-09, 4:2.5e-8, 5:4.5e-8, 6:5.79e-7}

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


if options.test in ('A','B','C','D'):  # set up sqrt geometry
   unique_xs=np.array(sorted(list(set(xCell[:]))))
   if min(unique_xs) >= 0.0:
      shiftx = unique_xs[3]  # put origin in a few cells to the right of the left edge of mesh
      xCell[:] = xCell[:] - shiftx
      xEdge[:] = xEdge[:] - shiftx
      xVertex[:] = xVertex[:] - shiftx


   # thickness
   x = 100.0e3 - np.absolute(xCell[:] - 100.0e3)
   thickness[0,:] = 6.0 *( np.sqrt( x[:] + 5.0e3) - np.sqrt(5.0e3) ) + 1.0
   thickness[0, x<0.0]=0.0  # make a margin

   # flat bed
   bedTopography[:] = 0.0

# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

if options.test == 'A':
  gridfile.variables['externalWaterInput'][:] = a_params[options.number] * 1000.0  # Convert from m/s to kg/m2/s

# melt
gridfile.variables['basalMeltInput'][:] = 0.0 * 1000.0 # no basal melting

# velocity
gridfile.variables['uReconstructX'][:] = 1.0e-6 # all tests use this velocity

# initial thickness
gridfile.variables['waterThickness'][0,:] = 0.05 * np.absolute(xCell[:] - 100.0e3)/100.0e3 # start with something to avoid weird adapative time steps initially

# initial waterPressure
gridfile.variables['waterPressure'][:] = 0.5 * 9.81 * 910.0 * thickness[:] * (0.97 + 0.03 * np.random.rand(thickness[:].size)) # start with half of Pice.  Last factor adds some random noise to disrupt symmetry.

gridfile.close()

print 'Successfully added initial conditions to: ', options.filename

