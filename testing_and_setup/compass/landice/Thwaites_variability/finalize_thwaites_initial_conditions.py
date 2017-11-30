#!/usr/bin/env python
'''
script to make final modifications needed for Thwaites initial condition.
run with one argument for the file to be modified.
'''

import sys
import netCDF4
import numpy as np

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file in which to finalize Thwaites initial conditions", metavar="FILE")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'

f = netCDF4.Dataset(options.filename,'r+')
thickness = f.variables['thickness'][0,:]
sfcMassBal = f.variables['sfcMassBal'][0,:]
bedTopography = f.variables['bedTopography'][0,:]
beta = f.variables['beta'][0,:]
cOnC= f.variables['cellsOnCell'][:]
ncOnC= f.variables['nEdgesOnCell'][:]
nCells = len(f.dimensions['nCells'])
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]


# Remove one-cell embayments from the calving front.
# These can lead to unrealistically large velocities for some reason.
# (better solution would be to understand why that is and stop it.)
for iter in range(3):  # iterate 3x to make sure all such locations have been filled in
   thktemp = np.copy(thickness)
   for c in range(nCells):
      if thickness[c] == 0.0:
         oceanNeighb = 0
         for n in range(ncOnC[c]):
            if thickness[cOnC[c,n]-1] == 0.0:
               oceanNeighb += 1
         if oceanNeighb == 1:
            thktemp[c] = thickness[cOnC[c,:ncOnC[c]]-1].mean()
   thickness = np.copy(thktemp)
f.variables['thickness'][0,:] = thickness
# put in large negative SMB beyond current calving front
killVal = -10.0
sfcMassBal[thickness<1.0e-3] = killVal
f.variables['sfcMassBal'][0,:] = sfcMassBal


# put small-ish beta value beyond current GL in case GL advances
beta[thickness*910.0/1028.0+bedTopography < 0.0] = 200.0

# optional: put larger beta on  shelf pinning point
#beta[((yCell - (-453466.0))**2 + (xCell - (-1.5925e6))**2)**0.5 < 12000.0] = 250.0

f.variables['beta'][0,:] = beta

f.close()
