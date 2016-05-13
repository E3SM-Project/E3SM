#!/usr/bin/env python
# This script sets up MISMIP3D Stnd experiment.
# see http://homepages.ulb.ac.be/~fpattyn/mismip3d/Mismip3Dv12.pdf

import sys
from netCDF4 import Dataset
from math import sqrt
import numpy as np
from collections import Counter

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file in which to setup MISMIP3D Stnd", metavar="FILE")
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
    if nVertLevels != 10:
         print 'nVerLevels in the supplied file was ', nVertLevels, '.  10 levels is a preliminary value to be used  with this test case.'
except:
    sys.exit('Error: The grid file specified is missing needed dimensions.')



# put the domain origin in the center of the center cell in the y-direction and in the 2nd row on the x-direction
# Only do this if it appears this has not already been done:
xCell = gridfile.variables['xCell'][:]
yCell = gridfile.variables['yCell'][:]
if yCell.min() > 0.0:
   print 'Shifting domain origin, because it appears that this has not yet been done.'
   unique_ys=np.array(sorted(list(set(yCell[:]))))
   targety = (unique_ys.max() - unique_ys.min()) / 2.0 + unique_ys.min()  # center of domain range
   best_y=unique_ys[ np.absolute((unique_ys - targety)) == np.min(np.absolute(unique_ys - (targety))) ][0]
   print 'Found a best y value to use of:' + str(best_y)
   
   unique_xs=np.array(sorted(list(set(xCell[:]))))
   targetx = (unique_xs.max() - unique_xs.min()) / 2.0 + unique_xs.min()  # center of domain range
   best_x=unique_xs[ np.absolute((unique_xs - targetx)) == np.min(np.absolute(unique_xs - (targetx))) ][0]
   print 'Found a best x value to use of:' + str(best_x)
   
   xShift = -1.0 * best_x
   yShift = -1.0 * best_y
   gridfile.variables['xCell'][:] = xCell + xShift
   gridfile.variables['yCell'][:] = yCell + yShift
   xCell = xCell + xShift
   yCell = yCell + yShift
   gridfile.variables['xEdge'][:] = gridfile.variables['xEdge'][:] + xShift
   gridfile.variables['yEdge'][:] = gridfile.variables['yEdge'][:] + yShift
   gridfile.variables['xVertex'][:] = gridfile.variables['xVertex'][:] + xShift
   gridfile.variables['yVertex'][:] = gridfile.variables['yVertex'][:] + yShift

   # Need to adjust geometry along top and bottom boundaries to get flux correct there
   # Essentially, we only want to model the interior half of those cells
   # Adding this here because we only want to do this if it hasn't been done before.
   print "Adjusting areaCell and dvEdge for cells along north and south boundaries"

   # This method is assuming a periodic_hex mesh!

   # Adjust area in half for N/S boundary cells
   unique_ys=np.array(sorted(list(set(yCell[:]))))  # recalculate after above adjustment
   areaCell = gridfile.variables['areaCell']  # Note just getting object here
   areaCell[ np.nonzero(yCell == unique_ys[0]) ] *= 0.5  # cut area in half for south row
   areaCell[ np.nonzero(yCell == unique_ys[-1]) ] *= 0.5  # cut area in half for north row

   # Adjust length in half for edges connecting N/S boundary cells
   yEdge = gridfile.variables['yEdge'][:]
   dvEdge = gridfile.variables['dvEdge']  # Note just getting object here
   unique_ys_edge = np.array(sorted(list(set(yEdge[:]))))
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[0]) ] *= 0.0  # zero out the edges on the boundary (not necessary cause velocity will also be zero)
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[1]) ] *= 0.5  # cut length in half for edges between boundary cells
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[-1]) ] *= 0.0  # zero out the edges on the boundary (not necessary cause velocity will also be zero)
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[-2]) ] *= 0.5  # cut length in half for edges between boundary cells


#   print np.array(sorted(list(set(yCell[:]))))

# bed slope defined by b(m)=-100km-x(km)
print "Defining bed topography."
topg = np.zeros((1,nCells,))
topg[0,np.nonzero(xCell>=0.0)]= -100.0 - xCell[np.nonzero(xCell>=0.0)]/1000.0
topg[0,np.nonzero(xCell<0.0)] = -100.0 + xCell[np.nonzero(xCell< 0.0)]/1000.0
gridfile.variables['bedTopography'][:] = topg[:]

# SMB
print "Defining SMB."
SMB = np.zeros((1,nCells,))
# 0.5 m/yr is the standard value.  0.3 m/yr is also tested in MISMIP3D.
# Convert from units of m/yr to kg/m2/s using appropriate ice density
SMB[:] = 0.5 *900.0/(3600.0*24.0*365.0)
# Add a 'gutter' along the eastern edge
SMB[ 0,np.nonzero(xCell >  800000.0) ] = -100.0
SMB[ 0,np.nonzero(xCell < -800000.0) ] = -100.0
gridfile.variables['sfcMassBal'][:] = SMB[:]

# Thickness initial condition is no ice.
print "Defining thickness."
thickness = np.zeros((nCells,))

thicknessICtype = 1  # 1=uniform, 2=b.l. solution
if thicknessICtype == 1:
   thickness[:] = 200.0 # 10m is the IC prescribed by the Stnd experiment, but nothing is grounded at that thickness, so make it thick enough to ground
   thickness[ np.nonzero(np.absolute(xCell) > 800000.0) ] = 0.0
   gridfile.variables['thickness'][0,:] = thickness[:]
elif thicknessICtype == 2:
   # Integrate Eq. 25 from Schoof 2007 (JGR) from a chosen xg position
   # Need to choose G.L. position.  Best to not choose the actual solution because at coarse resolution the GL won't move much and we'll get the right answer only because we started there...
   #xg = 500000.1  # m  This is approximately where it should end up, but better not to use this i.c. for the actual experiment.  This is primarily good for testing it doesn't move away from here...
   xg = 400000.1  # m
   C = 1.0e7    # Pa m^-1/3 s^1/3
   m = 1.0/3.0
   a = 0.3/3.14e7  # m/yr converted to m/s
   rhoi = 900.0  # kg/m3
   rhow = 1000.0 # kg/m3
   g = 9.81  #m/s2
   # calculations:
   ind = np.logical_and(xCell>=0.0, xCell<=xg)  # indices where xCell is between divide and GL
   unique_xs=np.array(sorted(list(set(xCell[ind])), reverse=True))  # returns a list of x values from GL to divide (descending)
   print unique_xs
#   i0 = np.nonzero(xCell == 0.0)[0][0] # index at divide
   hg = rhow/rhoi * -1.0*(-100.0 - xg/1000.0)  # thickness at GL
   print 'xg, hg=', xg, hg

   ig = np.nonzero(np.logical_or(xCell >= xg, xCell <= -xg))[0]  # indices at GL cells and shelf cells
   thickness[ig] = hg/2.0   # make shelf thickness same as GL halved - this is reasonably close to the shelf thickness at SS and avoids large velocities (and small CFL time step limits if using the full GL thickness)
   previous_x = xg
   previous_h = hg

   dbdx = 1.0/1000.0  # bed slope is uniform
   for x in unique_xs:
      q = a * x
      # first order integration
      h = previous_h + (dbdx - C/(rhoi*g) * abs(q)**(m-1.0) * q / previous_h**(m+1.0)) * (x - previous_x)
      if h != h:
          h = previous_h  # Nan at divide
      these_x = np.nonzero(xCell == x)[0]  # get the indices we are calculating for
      other_x = np.nonzero(xCell == -x)[0]  # get the mirror side indices
#      print x, these_x, xCell[these_x], h
      thickness[these_x] = h
      thickness[other_x] = h
      previous_x = x
      previous_h = h
   thickness[thickness != thickness] = 50000.0  # NaN check!
   thickness[ np.nonzero(np.absolute(xCell) > 800000.0) ] = 0.0
   gridfile.variables['thickness'][0,:] = thickness[:]



# For now approximate boundary conditions with 0 velocity.
print "Defining velocity boundary conditions."
# This is not correct.
# west boundary should be dh/dx=ds/dx=0.
# north and south boundaries should be no slip lateral boundaries.
# Dirichlet velocity mask
kinbcmask = np.zeros((1, nCells, nVertInterfaces))
kinbcmask[:, np.nonzero(yCell == yCell.min()), : ] = 1 # south row
kinbcmask[:, np.nonzero(yCell == yCell.max()), : ] = 1 # north row
###kinbcmask[:, np.nonzero(xCell < 0.0), : ] = 1 # west boundary
gridfile.variables['dirichletVelocityMask'][:] = kinbcmask
# Dirichlet velocity values
gridfile.variables['uReconstructX'][:] = 0.0
gridfile.variables['uReconstructY'][:] = 0.0

# beta is not correct
print "Defining beta."
#gridfile.variables['beta'][:] = 1.0e7 / 3.14e7**(1.0/m)   # For the basal friction law being used, beta holds the 'C' coefficient in Pa m^-1/3 s^1/3
gridfile.variables['beta'][:] = 31880.0  # For the basal friction law being used, beta holds the 'C' coefficient.  The beta units in MPAS are a mess right now.  This value translates to 10^7 Pa m^-1/3 s^1/3

# Setup layerThicknessFractions
gridfile.variables['layerThicknessFractions'][:] = 1.0 / float(nVertLevels)

gridfile.sync()
gridfile.close()

print 'Successfully added MISMIP3D initial conditions to: ', options.filename

