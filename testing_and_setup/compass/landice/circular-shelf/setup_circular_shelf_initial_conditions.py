#!/usr/bin/env python
# This script runs a "Circular Shelf Experiment".

import sys, numpy
from netCDF4 import Dataset
from math import sqrt

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file in which to setup circular shelf", metavar="FILE")
parser.add_option("-b","--beta", dest="use_beta", action="store_true", help="Use this flag to use a high value of the field 'beta' to specify 'no-slip' conditions in the grounded portion of the domain.  The default is to use Dirichlet boundary conditions.")
parser.add_option("-7","--7cells", dest="use_7cells", action="store_true", help="Use this flag to create 7 grounded cells (the center cell and its 6 neigbors).  The default is to only ground the center cell.")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'


# Open the file, get needed dimensions
gridfile = Dataset(options.filename,'r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
if nVertLevels != 5:
     print 'nVertLevels in the supplied file was ', nVertLevels, '.  This test case is typically run with 5 levels.'
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
cellsOnCell= gridfile.variables['cellsOnCell']
# Get b.c. variables
SMB = gridfile.variables['sfcMassBal']



# Center the dome in the center of the cell that is closest to the center of the domain.
# Only do this if it appears this has not already been done:
if xVertex[:].min() == 0.0:
   print "Shifting x/y coordinates to center domain at 0,0."
   # Find center of domain
   x0 = xCell[:].min() + 0.5 * (xCell[:].max() - xCell[:].min() )
   y0 = yCell[:].min() + 0.5 * (yCell[:].max() - yCell[:].min() )
   # Calculate distance of each cell center from dome center
   r = ((xCell[:] - x0)**2 + (yCell[:] - y0)**2)**0.5
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
# Now update our local values of the origin location and distance array (or assume these are correct because this grid has previously been shifted
x0 = 0.0
y0 = 0.0
r = ((xCell[:] - x0)**2 + (yCell[:] - y0)**2)**0.5
centerCellIndex = numpy.abs(r[:]).argmin()

# Make a circular ice mass
# Define dome dimensions - all in meters
r0 = 21000.0
thickness[:] = 0.0  # initialize to 0.0
# Calculate the dome thickness for cells within the desired radius (thickness will be NaN otherwise)
thickness_field = thickness[0,:]
thickness_field[r<r0] = 1000.0
thickness[0,:] = thickness_field

# flat bed at -2000 m everywhere with a single grounded point
bedTopography[:] = -2000.0  
bedTopography[0, centerCellIndex] = -880.0
if options.use_7cells:
   print 'Making the grounded portion of the domain cover 7 cells - the center cell and its 6 neighbors.'
   bedTopography[0, cellsOnCell[centerCellIndex,:]-1] = -880.0  # use this to make the grounded area 7 cells instead of 1
else:
   print 'Making the grounded portion of the domain cover 1 cell - the center cell.'

if options.use_beta:
   print 'Setting no-slip on the grounded portion of the domain by setting a high beta field there.'
   beta = gridfile.variables['beta']
   # beta is 0 everywhere except a high value in the grounded cell
   beta[:] = 0.
   beta[centerCellIndex] = 1.0e8
   if options.use_7cells:
      beta[cellsOnCell[centerCellIndex,:]-1] = 1.0e8 # use this to make the grounded area 7 cells instead of 1
else: # use Dirichlet b.c.
   print 'Setting no-slip on the grounded portion of the domain by setting no-slip Dirichlet velocity boundary conditions there.'
   dirMask = gridfile.variables['dirichletVelocityMask']
   uvel = gridfile.variables['uReconstructX']
   vvel = gridfile.variables['uReconstructY']
   dirMask[:] = 0
   uvel[:] = 0.0
   vvel[:] = 0.0
   # Apply mask to basal level of grounded portion only
   dirMask[:, centerCellIndex, -1] = 1
   if options.use_7cells:
      dirMask[:, cellsOnCell[centerCellIndex,:]-1, -1] = 1

# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels
# boundary conditions
SMB[:] = 0.0  # m/yr
# Convert from units of m/yr to kg/m2/s using an assumed ice density
SMB[:] = SMB[:] *910.0/(3600.0*24.0*365.0)

gridfile.close()

print '\nSuccessfully added circular-shelf initial conditions to: ', options.filename



