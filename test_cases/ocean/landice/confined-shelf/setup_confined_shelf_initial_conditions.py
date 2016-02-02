#!/usr/bin/env python
# This script sets up a "Confined Shelf Experiment".
# see http://homepages.vub.ac.be/~phuybrec/eismint/shelf-descr.pdf

import sys
from netCDF4 import Dataset
from math import sqrt
import numpy as np
from collections import Counter

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file in which to setup confined shelf", metavar="FILE")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'


# Open the file, get needed dimensions
try:
    gridfile = Dataset(options.filename,'r+')
    nCells = len(gridfile.dimensions['nCells'])
    nVertLevels = len(gridfile.dimensions['nVertLevels'])
    maxEdges = len(gridfile.dimensions['maxEdges'])
    if nVertLevels != 5:
         print 'nVerLevels in the supplied file was ', nVertLevels, '.  5 levels are typically used with this test case.'
    # Get variables
    xCell = gridfile.variables['xCell']
    yCell = gridfile.variables['yCell']
    xEdge = gridfile.variables['xEdge']
    yEdge = gridfile.variables['yEdge']
    xVertex = gridfile.variables['xVertex']
    yVertex = gridfile.variables['yVertex']
    thickness = gridfile.variables['thickness']
    bedTopography = gridfile.variables['bedTopography']
    beta = gridfile.variables['beta']
    layerThicknessFractions = gridfile.variables['layerThicknessFractions']
    temperature = gridfile.variables['temperature']
    cellsOnCell= gridfile.variables['cellsOnCell']
    # Get b.c. variables
    kinbcmask = gridfile.variables['dirichletVelocityMask']
    uvel = gridfile.variables['uReconstructX']
    vvel = gridfile.variables['uReconstructY']
    SMB = gridfile.variables['sfcMassBal']
except:
    sys.exit('Error: The grid file specified is either missing or lacking needed dimensions/variables.')



# put the domain origin in the center of the center cell in the x-direction and in the 2nd row on the y-direction
# Only do this if it appears this has not already been done:
if xVertex[:].min() == 0.0:
   print 'Shifting domain origin to center of shelf front, because it appears that this has not yet been done.'
   unique_xs=np.array(sorted(list(set(xCell[:]))))
   targetx = (unique_xs.max() - unique_xs.min()) / 2.0 + unique_xs.min()  # center of domain range
   best_x=unique_xs[ np.absolute((unique_xs - targetx)) == np.min(np.absolute(unique_xs - (targetx))) ][0]
   print 'Found a best x value to use of:' + str(best_x)
   
   unique_ys=np.array(sorted(list(set(yCell[:]))))
#   print unique_ys
   best_y = unique_ys[5]  # get 6th value
   print 'Found a best y value to use of:' + str(best_y)
   
   xShift = -1.0 * best_x
   yShift = -1.0 * best_y
   xCell[:] = xCell[:] + xShift
   yCell[:] = yCell[:] + yShift
   xEdge[:] = xEdge[:] + xShift
   yEdge[:] = yEdge[:] + yShift
   xVertex[:] = xVertex[:] + xShift
   yVertex[:] = yVertex[:] + yShift

#   print np.array(sorted(list(set(yCell[:]))))

# Make a square ice mass
# Define square dimensions - all in meters
L = 200000.0

shelfMask = np.logical_and( 
                  np.logical_and(xCell[:]>=-L/2.0, xCell[:]<=L/2.0), 
                  np.logical_and(yCell[:]>=0.0, yCell[:]<=L) ) 
# now grow it by one cell
shelfMaskWithGround = np.zeros( (nCells,), dtype=np.int16)
for c in range(nCells):
    if shelfMask[c] == 1:
        for n in range(maxEdges):
            neighbor = cellsOnCell[c,n] - 1 # fortran to python indexing
            if neighbor >= 0:
                shelfMaskWithGround[neighbor] = 1
# but remove the south side extension
shelfMaskWithGround[ np.nonzero(yCell[:]<0.0)[0] ] = 0

thickness[:] = 0.0  # initialize to 0.0
thickness[0, np.nonzero(shelfMaskWithGround==1)[0] ] = 500.0

# flat bed at -2000 m everywhere but grounded around the edges
bedTopography[np.nonzero(shelfMask==1)[0]] = -2000.0
bedTopography[np.nonzero(shelfMask==0)[0]] = -440.0

# Dirichlet velocity mask
kinbcmask[:] = 0
#kinbcmask[:, np.nonzero(shelfMask==0)[0], :] = 1
kinbcmask[:, np.nonzero(np.logical_and(shelfMask==0, shelfMaskWithGround==1))[0], :] = 1
#kinbcmask[:, np.nonzero(yCell[:]<0.0)[0] ] = 0
# Need to extend this mask south by one cell so that the extended FEM mask will still have the 0 velo on the edges...
theSides = Counter(xCell[ np.nonzero(kinbcmask[0,:])[0] ]).most_common(4)  # need the 4 most common x positions.
for side in theSides:
    thesideindices = np.nonzero( np.logical_and( xCell[:] == side[0] , yCell[:] <= 0.0 ) )[0]
    kinbcmask[:, thesideindices] = 1

# Dirichlet velocity values
uvel[:] = 0
vvel[:] = 0

# beta is 0 everywhere (strictly speaking it should not be necessary to set this)
beta[:] = 0.

# constant, arbitrary temperature, degrees C
temperature[:] = 273.15 

# Setup layerThicknessFractions
layerThicknessFractions[:] = 1.0 / nVertLevels

# boundary conditions
SMB[:] = 0.0  # m/yr
# Convert from units of m/yr to kg/m2/s using an assumed ice density
SMB[:] = SMB[:] *910.0/(3600.0*24.0*365.0)

gridfile.close()

print 'Successfully added confined-shelf initial conditions to: ', options.filename



