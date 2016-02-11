#!/usr/bin/env python
# Generate initial conditions for EISMINT-1 moving margin land ice test case

import sys, numpy
from netCDF4 import Dataset as NetCDFFile
from math import sqrt
import numpy

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file to setup dome", metavar="FILE")
options, args = parser.parse_args()

if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'


# Open the file, get needed dimensions
try:
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
    temperature = gridfile.variables['temperature']
    # Get b.c. variables
    SMB = gridfile.variables['sfcMassBal']
except:
    sys.exit('Error: The grid file specified is either missing or lacking needed dimensions/variables.')

# Find center of domain
x0 = xCell[:].min() + 0.5 * (xCell[:].max() - xCell[:].min() )
y0 = yCell[:].min() + 0.5 * (yCell[:].max() - yCell[:].min() )
# Calculate distance of each cell center from dome center
r = ((xCell[:] - x0)**2 + (yCell[:] - y0)**2)**0.5

# Center the dome in the center of the cell that is closest to the center of the domain.
#   NOTE: for some meshes, maybe we don't want to do this - could add command-line argument controlling this later.
putOriginOnACell = True
if putOriginOnACell:
   centerCellIndex = numpy.abs(r[:]).argmin()
   #print x0, y0, centerCellIndex, xCell[centerCellIndex], yCell[centerCellIndex]
   xShift = -1.0 * xCell[centerCellIndex]
   yShift = -1.0 * yCell[centerCellIndex]
   xCell[:] = xCell[:] + xShift
   yCell[:] = yCell[:] + yShift
   xEdge[:] = xEdge[:] + xShift
   yEdge[:] = yEdge[:] + yShift
   xVertex[:] = xVertex[:] + xShift
   yVertex[:] = yVertex[:] + yShift
   # Now update origin location and distance array
   x0 = 0.0
   y0 = 0.0
   r = ((xCell[:] - x0)**2 + (yCell[:] - y0)**2)**0.5


# Assign variable values for EISMINT-1 experiment
# Start with no ice
thickness[:] = 0.0
# flat bed at sea level
bedTopography[:] = 0.0
# constant, arbitrary temperature, degrees C
temperature[:] = 0.0 
# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

# boundary conditions
# SMB field specified by EISMINT, constant in time for EISMINT-1
# Define function for SMB
rhoi = 910.0
scyr = 3600.0*24.0*365.0
Rel = 450.0  # km
s = 10.0**-2     # given in units of m/a/km, 
smb=numpy.minimum(0.5, s * (Rel - r/1000.0)) # in m ice/yr
#s = s * rhoi / 1000.0 / (3600.0*24.0*365.0)  # converted to kg/m2/s/m using ice density of 910.0
#SMB[:] = numpy.minimum(0.5 * rhoi / (3600.0*24.0*365.0),   s * (Rel - d) )
SMB[0,:] = smb * rhoi / scyr  # in kg/m2/s

# Basal heat flux should be -42.e-3 once it is added.
# Surface temperature should be 270 K - 0.01 H when it is added.


gridfile.close()
print 'Successfully added initial conditions for EISMINT1-Moving Margin, experiment 1 to the file: ', options.filename

