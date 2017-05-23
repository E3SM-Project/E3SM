#!/usr/bin/env python
'''
Generate initial conditions for hydro-SHMIP land ice test case E
Details here: http://shmip.bitbucket.org/
'''

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import sys

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file to setup", metavar="FILE")
parser.add_option("-n", "--number", dest="number", type='int', help="test variant to set up, 1-5", metavar="NUMBER")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'

# Setup dictionaries of parameter values for each experiment
e_params = {1:0.05,  2:0.0, 3:-0.1, 4:-0.5, 5:-0.7}

# Open the file, get needed dimensions
gridfile = NetCDFFile(options.filename,'r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
nEdges = len(gridfile.dimensions['nEdges'])


# adjust mesh origin if it hasn't been done yet
yCell = gridfile.variables['yCell']
if yCell[:].min() >= 0.0:
   xCell = gridfile.variables['xCell']
   xEdge = gridfile.variables['xEdge']
   yEdge = gridfile.variables['yEdge']
   xVertex = gridfile.variables['xVertex']
   yVertex = gridfile.variables['yVertex']
   # x-origin should be a few cells in from the left side of domain
   xShift = np.unique(xEdge[:])[9]
   xCell[:] = xCell[:] - xShift
   xEdge[:] = xEdge[:] - xShift
   xVertex[:] = xVertex[:] - xShift

   # y origin should be in center of domain, but centered on a cell
   uniqueY = np.unique(yCell[:])
   yShift = uniqueY[len(uniqueY)//2]  # use middle cell center value
   yCell[:] = yCell[:] - yShift
   yEdge[:] = yEdge[:] - yShift
   yVertex[:] = yVertex[:] - yShift

   gridfile.sync()
   print "Shifted domain origin by (x,y):", xShift, yShift

xCell = gridfile.variables['xCell'][:]
yCell = gridfile.variables['yCell'][:]
xEdge = gridfile.variables['xEdge'][:]
yEdge = gridfile.variables['yEdge'][:]


# Geometry
para = e_params[options.number]
print "Setting up test number:", options.number
print "Using overdeepening parameter of:", para

# prescribed functions:
para_bench = 0.05
surface = lambda x: 100.0 * (x + 200.0)**0.25 + 1.0/60.0 * x - 2.0e10**0.25 + 1.0
f = lambda x,para: (surface(6.0e3) - para * 6.0e3) / 6.0e3**2 * x**2 + para * x
g = lambda y: 0.5e-6 * np.absolute(y)**3
h = lambda x, para: (-4.5 * x / 6.0e3 + 5.0) * (surface(x)-f(x, para)) / (surface(x) - f(x, para_bench) + np.finfo(float).eps)

# calculate geometry:
xClean = np.maximum(xCell, xCell*0.0)
#print np.unique(xClean)
upperSurface = surface(xClean)
#print upperSurface
bed = f(xClean, para) + g(yCell) * h(xClean, para)
thickness = upperSurface - bed
thickness[thickness<0.0]=0.0  # remove negative thickness
thickness[xCell<=0.0]=0.0 # remove the 1m min of ice in surface that remains when using xClean

# store values
gridfile.variables['bedTopography'][0,:] = bed
gridfile.variables['thickness'][0,:] = thickness
##gridfile.variables['sfcMassBal'][0,:] = upperSurface  # hijacking for debugging

gridfile.sync()


# Set up no flux mask along all edges
# loop over cells and check for ice on only one side
cONe = gridfile.variables['cellsOnEdge'][:]
waterFluxMask = gridfile.variables['waterFluxMask'][0,:]
waterFluxMask[:] = 0
dx = gridfile.variables['dcEdge'][:].min()
for i in range(nEdges):
   c1 = cONe[i,0]-1
   c2 = cONe[i,1]-1
   if int(thickness[c1] > 0.0) + int(thickness[c2] > 0.0) == 1: # only one neighboring cell has ice
      if xEdge[i]<dx and np.absolute(yEdge[i])<0.5*dx:   # leave an outflow region near terminus
         pass
      else:
         waterFluxMask[i] = 2
gridfile.variables['waterFluxMask'][0,:] = waterFluxMask
gridfile.sync()


# Setup layerThicknessFractions
gridfile.variables['layerThicknessFractions'][:] = 1.0 / nVertLevels

# Water input
waterInput = gridfile.variables['externalWaterInput'][0,:]
waterInput[:] = 0.0
waterInput[thickness[:]>0.0] = 1.158e-6 * 1000.0  # Convert from m/s to kg/m2/s
gridfile.variables['externalWaterInput'][0,:] = waterInput


# melt
gridfile.variables['basalMeltInput'][:] = 0.0 * 1000.0 # no basal melting

# velocity
gridfile.variables['uReconstructX'][:] = 1.0e-6 # all tests use this velocity

# initial water thickness
gridfile.variables['waterThickness'][0,:] = 0.005 * np.absolute(xCell[:] - 6.0e3)/6.0e3 # start with something very small but nonzero to avoid weird adapative time steps initially

# initial waterPressure
#gridfile.variables['waterPressure'][:] = 0.5 * 9.81 * 910.0 * thickness[:] * (0.97 + 0.03 * np.random.rand(thickness[:].size)) # start with half of Pice.  Last factor adds some random noise to disrupt symmetry.
gridfile.variables['waterPressure'][0,:] = 0.05 * 9.81 * 910.0 * thickness[:]  # start with something small to avoid unphysical gradients

gridfile.close()

print 'Successfully added initial conditions to: ', options.filename

