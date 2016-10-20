#!/usr/bin/env python
#WHL - Based on MISMIP3D setup script by Matt Hoffman, modified for MISMIP+.

# This script sets up initial conditions for the MISMIP+ experiment.
# See this paper for details:
#   X. Asay-Davis et al. (2015), Experimental design for three interrelated 
#   Marine Ice-Sheet and Ocean Model Intercomparison Projects, Geosci. Model Devel. Discuss.,
#   8, 9859-9924.
# Following grid setup and initialization, the code should be spun up for ~10-20 ka
#  without basal melting, which should result in a stable grounding line
#  crossing the center of the channel around x = 450 km.

import sys
from netCDF4 import Dataset
from math import sqrt
import numpy as np
from collections import Counter

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file in which to set up MISMIP+ Stnd", metavar="FILE")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'


# Open the file, get needed dimensions
try:
    gridfile = Dataset(options.filename,'r+')
    nCells = len(gridfile.dimensions['nCells'])
    nVertLevels = len(gridfile.dimensions['nVertLevels'])
    nVertInterfaces = nVertLevels + 1
    maxEdges = len(gridfile.dimensions['maxEdges'])
    #WHL- Maybe do production runs with fewer than 10 levels?
    if nVertLevels != 10:
         print 'nVertLevels in the supplied file was ', nVertLevels, '.  10 levels is a preliminary value to be used with this test case.'
except:
    sys.exit('Error: The grid file specified is missing needed dimensions.')


# Put the domain origin at the center of the lower left grid cell.
# Only do this if it appears this has not already been done:

#WHL - Need [:] to actually read the data into a numpy array, and read it just once. Good for performance.
xCell = gridfile.variables['xCell'][:]
yCell = gridfile.variables['yCell'][:]

if yCell.min() > 0.0:
   print 'Shifting domain origin, because it appears that this has not yet been done.'
   unique_xs = np.array(sorted(list(set(xCell[:]))))
   unique_ys = np.array(sorted(list(set(yCell[:]))))

   print 'unique_ys.min:', unique_ys.min()
   print 'unique_ys.max:', unique_ys.max()
   print 'unique_xs.min:', unique_xs.min()
   print 'unique_xs.max:', unique_xs.max()

   xShift = -1.0 * unique_xs.min()
   yShift = -1.0 * unique_ys.min()
   gridfile.variables['xCell'][:] = xCell + xShift
   gridfile.variables['yCell'][:] = yCell + yShift
   xCell = xCell + xShift
   yCell = yCell + yShift
   gridfile.variables['xEdge'][:] = gridfile.variables['xEdge'][:] + xShift
   gridfile.variables['yEdge'][:] = gridfile.variables['yEdge'][:] + yShift
   gridfile.variables['xVertex'][:] = gridfile.variables['xVertex'][:] + xShift
   gridfile.variables['yVertex'][:] = gridfile.variables['yVertex'][:] + yShift

   # Need to adjust geometry along top and bottom boundaries to get flux correct there.
   # Essentially, we only want to model the interior half of those cells.
   # Adding this here because we only want to do this if it hasn't been done before.
   # This method is assuming a periodic_hex mesh!
   print "Adjusting areaCell and dvEdge for cells along north and south boundaries"

   # Adjust area in half for N/S boundary cells
   unique_ys = np.array(sorted(list(set(yCell[:]))))  # recalculate after above adjustment
   areaCell = gridfile.variables['areaCell']  # Note just getting object here
   areaCell[ np.nonzero(yCell == unique_ys[0]) ] *= 0.5   # cut area in half for south row
   areaCell[ np.nonzero(yCell == unique_ys[-1]) ] *= 0.5  # cut area in half for north row

   # Adjust length in half for edges connecting N/S boundary cells
   yEdge = gridfile.variables['yEdge'][:]
   dvEdge = gridfile.variables['dvEdge']  # Note just getting object here
   unique_ys_edge = np.array(sorted(list(set(yEdge[:]))))
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[0]) ] *= 0.0  # zero out the edges on the boundary (not necessary because velocity will also be zero)
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[1]) ] *= 0.5   # cut length in half for edges between boundary cells
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[-1]) ] *= 0.0  # zero out the edges on the boundary (not necessary because velocity will also be zero)
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[-2]) ] *= 0.5  # cut length in half for edges between boundary cells


#The following function computes the MISMIP+ bed according to Asay-Davis et al. (2016)

def computeBed(x,y):
   x = x/1.e3     # m to km
   y = y/1.e3     # m to km
   B0 = -150.     # m
   B2 = -728.8    # m
   B4 = 343.91    # m
   B6 = -50.57    # m
   x_bar = 300.   # km
   x_tilde = x/x_bar
   dc = 500.      # m
   fc = 4.        # km
   wc = 24.       # km
   Ly = 80.       # km
   Bmax = -720.   # m
   B_x = B0 + B2*x_tilde**2 + B4*x_tilde**4 + B6*x_tilde**6
   B_y = dc / (1 + np.exp(-2*(y-Ly/2-wc)/fc)) + dc / (1 + np.exp(2*(y-Ly/2+wc)/fc))
   Bsum = B_x + B_y
   B = np.maximum(Bsum, Bmax)   # B >= Bmax
   return B
   
# Create the required variables in the netCDF file.

# Set bedTopography (this variable should always be present in the input file)
print "Defining bedTopography"

# Compute the bed topography
bedTopography = np.zeros((nCells,))
bedTopography = computeBed(xCell,yCell)
gridfile.variables['bedTopography'][0,:] = bedTopography[:]  # dimensions of gridfile variable are Time and nCells

# Debug: Print the topography along the bottom row.
#unique_xs = np.array(sorted(list(set(xCell[:]))))
#for iCell in range(1,nCells):
#   if yCell[iCell] == yCell.min():
#      if xCell[iCell] in unique_xs:
#         print xCell[iCell], bedTopography[iCell]

# Set the initial thickness
# Initial condition is uniform 100 m (except where x > xcalve)

print "Defining thickness"
xcalve = 640000.0       # m
init_thickness = 100.0  # m
thickness = np.zeros((nCells,))
for iCell in range(1,nCells):
   if xCell[iCell] < xcalve:
      thickness[iCell] = init_thickness

gridfile.variables['thickness'][0,:] = thickness[:]

# Set the surface mass balance.
# MISMIP+ assumes an SMB of 0.3 m/yr.
# Convert from m/yr to kg/m2/s using appropriate ice density.
# Assign a large negative SMB where x > xcalve, to prevent ice advancing.

print "Defining SMB"
SMB = np.zeros((nCells,))
rhoi = 918.0  # from Asay-Davis et al. (2016)
seconds_per_year = 3600.0 * 24.0 * 365.0
SMB[:] = 0.3 * rhoi/seconds_per_year
for iCell in range(1,nCells):
   if xCell[iCell] > xcalve:
      SMB[iCell] = -100.0

gridfile.variables['sfcMassBal'][0,:] = SMB[:]


# Approximate boundary conditions with a Dirichlet velocity mask (velocity = 0).
# Note: At the N and S boundaries, only the normal (y) velocity component
#       will be zeroed out in Albany.  The x component can be nonzero,
#       supporting a no-slip boundary condition.

if 'dirichletVelocityMask' in gridfile.variables:
   print 'dirichletVelocityMask already in gridfile'
   kinbcmask = gridfile.variables['dirichletVelocityMask']
else:
   print 'dirichletVelocityMask not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   dirichletVelocityMask = gridfile.createVariable('dirichletVelocityMask', datatype, ('Time','nCells','nVertInterfaces'))

print "Defining velocity boundary conditions"
kinbcmask = np.zeros((nCells, nVertInterfaces))
kinbcmask[np.nonzero(yCell == yCell.min()), : ] = 1 # south row
kinbcmask[np.nonzero(yCell == yCell.max()), : ] = 1 # north row
kinbcmask[np.nonzero(xCell < 0.0), : ] = 1          # west boundary
gridfile.variables['dirichletVelocityMask'][0,:] = kinbcmask

# Set the initial velocities to zero to enforce Dirichlet BC..
# May not be necessary, but doing this to be on the safe side.
if 'uReconstructX' in gridfile.variables:
   print 'uReconstructX already in gridfile'
else:
   print 'uReconstructX not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   uReconstructX = gridfile.createVariable('uReconstructX', datatype, ('Time','nCells'))

if 'uReconstructY' in gridfile.variables:
   print 'uReconstructY already in gridfile'
else:
   print 'uReconstructY not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   uReconstructY = gridfile.createVariable('uReconstructY', datatype, ('Time','nCells'))

gridfile.variables['uReconstructX'][0,:] = 0.0
gridfile.variables['uReconstructY'][0,:] = 0.0


# Set basal traction coefficient, beta.
# For now, assume a Weertman-type power law, tau_b = C * u^(1/m), where C = beta.
# Asay-Davis et al. (2016) specify C = 3.160 x 10^6 Pa m^{-1/3} s^{1/3} for power-law friction,
#  with friction-law exponent m = 3.
# Later, we could support a Tsai friction law.

if 'beta' in gridfile.variables:
   print 'beta already in gridfile'
   kinbcmask = gridfile.variables['beta']
else:
   print 'beta not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   beta = gridfile.createVariable('beta', datatype, ('Time','nCells'))

print "Defining beta"
# For the Weertman power law, beta holds the 'C' coefficient.  The beta units in MPAS are a mess right now.  
# In the MISMIP3D setup script, C = 10^7 Pa m^-1/3 s^1/3 translates to beta = 31880.
# For MISMIP+, C = 3.160 x 10^6 Pa m^-1/3 s^1/3 translates to beta = 10002.
 
C = 3.160e6   # Pa m^{-1/3} s^{1/3}
C = C / seconds_per_year**(1.0/3.0)  # convert to MPAS units
gridfile.variables['beta'][0,:] = C


# Set up layerThicknessFractions
gridfile.variables['layerThicknessFractions'][:] = 1.0 / float(nVertLevels)

#WHL - sync ensures that values are written to the file before closing.
gridfile.sync()
gridfile.close()

print 'Successfully added MISMIP+ initial conditions to: ', options.filename
