#!/usr/bin/env python
# Generate initial conditions for EISMINT-2 A land ice test case
# Test case is described in:
# Payne, A. J., Huybrechts, P., Calov, R., Fastook, J. L., Greve, R., Marshall, S. J., Marsiat, I., Ritz, C., Tarasov, L. and Thomassen, M. P. A.: Results from the EISMINT model intercomparison: the effects of thermomechanical coupling, J. Glaciol., 46(153), 227-238, 2000.

import sys, os
from netCDF4 import Dataset as NetCDFFile
from math import sqrt
import numpy
import shutil, glob

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-e", "--experiment", dest="exp", type='string', help="experiment to setup", metavar="EXP")
options, args = parser.parse_args()

if not options.exp:
    sys.exit('Error: No experiment specified.  Please specify an experiment to setup with the -e option')
experiment = options.exp.lower()  # allow lower or upper case to be input, but use lower case in the script logic
if experiment in ('a','b','c','d','f','g'):
    print 'Setting up EISMINT2 Experiment ' + experiment
else:
    sys.exit("Error: Invalid experiment specified.  Please specify an experiment between 'a' and 'g', excluding 'e'")


# Setup dictionaries of parameter values for each experiment
# Mmax: Maximum SMB at center of domain (m a-1)
# Sb: gradient of SMB with horizontal distance (m a-1 km-1)
# Rel: radial distance from summit where SMB = 0 (km)
# Tmin: surface temperature at summit (K)
# ST: gradient of air temperature with horizontal distance (K km-1)
# beta: basal traction coefficient (Pa m-1 a)
#       Note: beta is the inverse of parameter B in Payne et al. (2000)
a_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':450.0, 'Tmin':238.15, 'ST':1.67e-2, 'beta':1.0e8}
b_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':450.0, 'Tmin':243.15, 'ST':1.67e-2, 'beta':1.0e8}
c_params = {'Mmax':0.25, 'Sb':10.0**-2, 'Rel':425.0, 'Tmin':238.15, 'ST':1.67e-2, 'beta':1.0e8}
d_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':425.0, 'Tmin':238.15, 'ST':1.67e-2, 'beta':1.0e8}
f_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':450.0, 'Tmin':223.15, 'ST':1.67e-2, 'beta':1.0e8}
g_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':450.0, 'Tmin':238.15, 'ST':1.67e-2, 'beta':1.0e3}
xsummit = 750000.0; ysummit = 750000.0
rhoi = 910.0
scyr = 3600.0*24.0*365.0

# Setup dictionary of dictionaries for each experiment
exp_params = {'a':a_params, 'b':b_params, 'c':c_params, 'd':d_params, 'f':f_params, 'g':g_params}

# Some experiments start from scratch, others start from the SS of a previous experiment
#filename = 'eismint2' + experiment + '.input.nc'
filename = 'landice_grid.nc'
# experiments that start from scratch
if experiment in ('a', 'f', 'g'):
    # we will build the mesh from scratch
    shutil.copyfile('../setup_mesh/landice_grid.nc', filename)
else:
    # use the final state of experiment A
    try:
      stat = os.system('ncks -O -d Time,-1 ../experiment_A/output.nc ./' + filename)
      if stat != 0:
         raise error('ncks error')
    except:
      sys.exit("Error: problem building initial condition file from final state of Experiment A.  Make sure Experiment A has completed successfully before running tests B, C, or D.")

# Open the new input file, get needed dimensions & variables
gridfile = NetCDFFile(filename,'r+')
nVertLevels = len(gridfile.dimensions['nVertLevels'])
# Get variables
xCell = gridfile.variables['xCell'][:]
yCell = gridfile.variables['yCell'][:]
xEdge = gridfile.variables['xEdge'][:]
yEdge = gridfile.variables['yEdge'][:]
xVertex = gridfile.variables['xVertex'][:]
yVertex = gridfile.variables['yVertex'][:]

# ===================
# initial conditions
# ===================
# If starting from scratch, setup dimension variables and initial condition variables
if experiment in ('a', 'f', 'g'):
    # Find center of domain
    x0 = xCell[:].min() + 0.5 * (xCell[:].max() - xCell[:].min() )
    y0 = yCell[:].min() + 0.5 * (yCell[:].max() - yCell[:].min() )
    # Calculate distance of each cell center from dome center
    r = ((xCell[:] - x0)**2 + (yCell[:] - y0)**2)**0.5

    # Center the dome in the center of the cell that is closest to the center of the domain.
    centerCellIndex = numpy.abs(r[:]).argmin()
    # EISMINT-2 puts the center of the domain at 750,750 km instead of 0,0.
    # Adjust to use that origin.
    #print x0, y0, centerCellIndex, xCell[centerCellIndex], yCell[centerCellIndex]
    xShift = -1.0 * xCell[centerCellIndex] + xsummit
    yShift = -1.0 * yCell[centerCellIndex] + ysummit
    xCell[:] = xCell[:] + xShift
    yCell[:] = yCell[:] + yShift
    xEdge[:] = xEdge[:] + xShift
    yEdge[:] = yEdge[:] + yShift
    xVertex[:] = xVertex[:] + xShift
    yVertex[:] = yVertex[:] + yShift
    gridfile.variables['xCell'][:] = xCell[:]
    gridfile.variables['yCell'][:] = yCell[:]
    gridfile.variables['xEdge'][:] = xEdge[:]
    gridfile.variables['yEdge'][:] = yEdge[:]
    gridfile.variables['xVertex'][:] = xVertex[:]
    gridfile.variables['yVertex'][:] = yVertex[:]

    # Assign initial condition variable values for EISMINT-2 experiment
    # Start with no ice
    gridfile.variables['thickness'][:] = 0.0
    # flat bed at sea level
    gridfile.variables['bedTopography'][:] = 0.0
    # constant, arbitrary temperature, degrees K (doesn't matter since there is no ice initially)
    gridfile.variables['temperature'][:] = 273.15
    # Setup layerThicknessFractions
    gridfile.variables['layerThicknessFractions'][:] = 1.0 / nVertLevels
else:
    StrLen= len(gridfile.dimensions['StrLen'])
    gridfile.variables['xtime'][0,:] = list('000000-01-01_00:00:00'.ljust(StrLen, ' '))

# Now update/set origin location and distance array
r = ((xCell[:] - xsummit)**2 + (yCell[:] - ysummit)**2)**0.5

# ===================
# boundary conditions
# ===================
# Define values prescribed by Payne et al. 2000 paper.

params = exp_params[experiment]
print "Parameters for this experiment:", params

# SMB field specified by EISMINT, constant in time for EISMINT2
# It is a function of geographical position (not elevation)
Mmax = params['Mmax'] / scyr  # maximum accumulation rate [m/yr] converted to [m/s]
Sb = params['Sb'] / scyr / 1000.0  # gradient of accumulation rate change with horizontal distance  [m/a/km] converted to [m/s/m]
Rel = params['Rel'] * 1000.0  # accumulation rate at 0 position  [km] converted to [m]

SMB = numpy.minimum(Mmax, Sb * (Rel - r)) # [m ice/s]
SMB = SMB * rhoi  # in kg/m2/s
if 'sfcMassBal' in gridfile.variables:
   sfcMassBalVar = gridfile.variables['sfcMassBal']
else:
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   sfcMassBalVar = gridfile.createVariable('sfcMassBal', datatype, ('Time', 'nCells'))
sfcMassBalVar[0,:] = SMB


# Surface temperature
Tmin = params['Tmin']  # minimum surface air temperature [K]
ST = params['ST'] / 1000.0  # gradient of air temperature change with horizontal distance [K/km] converted to [K/m]
if 'surfaceAirTemperature' in gridfile.variables:
   surfaceAirTemperatureVar = gridfile.variables['surfaceAirTemperature']
else:
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   surfaceAirTemperatureVar = gridfile.createVariable('surfaceAirTemperature', datatype, ('Time', 'nCells'))
surfaceAirTemperatureVar[0,:] = Tmin + ST * r

# beta
beta = params['beta']
if 'beta' in gridfile.variables:
   betaVar = gridfile.variables['beta']
else:
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   betaVar = gridfile.createVariable('beta', datatype, ('Time', 'nCells'))
betaVar[0,:] = beta

gridfile.close()
print 'Successfully added initial conditions for EISMINT2, experiment '+experiment+' to the file: ', filename

