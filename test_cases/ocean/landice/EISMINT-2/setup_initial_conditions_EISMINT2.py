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
if experiment in ('a','b','c','d','f'):
    print 'Setting up EISMINT2 Experiment ' + experiment
else:
    sys.exit("Error: Invalid experiment specified.  Please specify an experiment between 'a' and 'f', excluding 'd'")

# Move to requested directory
try:
   os.chdir("experiment_" + experiment)
except:
   sys.exit("Error: unable to move to experiment directory experiment_" + experiment)



# Setup dictionaries of parameter values for each experiment
a_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':450.0, 'Tmin':238.15, 'ST':1.67e-2}
b_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':450.0, 'Tmin':243.15, 'ST':1.67e-2}
c_params = {'Mmax':0.25, 'Sb':10.0**-2, 'Rel':425.0, 'Tmin':238.15, 'ST':1.67e-2}
d_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':425.0, 'Tmin':238.15, 'ST':1.67e-2}
f_params = {'Mmax':0.5,  'Sb':10.0**-2, 'Rel':450.0, 'Tmin':223.15, 'ST':1.67e-2}
xsummit = 750000.0; ysummit = 750000.0
rhoi = 910.0
scyr = 3600.0*24.0*365.0

# Setup dictionary of dictionaries for each experiment
exp_params = {'a':a_params, 'b':b_params, 'c':c_params, 'd':d_params, 'f':f_params}

# Some experiments start from scratch, others start from the SS of a previous experiment
filename = 'eismint2' + experiment + '.input.nc'
# experiments that start from scratch
if experiment in ('a', 'f'):
    # get a empty land ice grid to use from the parent directory
    try:
      shutil.copy("../landice_grid.nc", filename)
      for file in glob.glob('../graph.info*'):
          shutil.copy(file, '.')
    except:
      sys.exit("Error: problem copying ../landice_grid.nc and/or graph.info* to this directory")
else:
    # use the final state of experiment A
    try:
      stat = os.system('ncks -O -d Time,-1 ../experiment_a/eismint2a.output.nc ./' + filename) 
      if stat != 0:
         raise error('ncks error')
    except:
      sys.exit("Error: problem building initial condition file from final state of Experiment A.  Make sure Experiment A has completed successfully before running tests B, C, or D.")

# setup graph.info files
try:
  for filepath in glob.glob('../graph.info.part.*'):
      #os.symlink(filepath, os.path.basename(filepath))
      shutil.copy(filepath, '.')
except:
  sys.exit("Error: problem copying graph.info.part.* files to experiment directory")

# Open the new input file, get needed dimensions & variables
try:
    gridfile = NetCDFFile(filename,'r+')
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
    normalVelocity = gridfile.variables['normalVelocity']
    layerThicknessFractions = gridfile.variables['layerThicknessFractions']
    temperature = gridfile.variables['temperature']
    # Get b.c. variables
    SMB = gridfile.variables['sfcMassBal']
except:
    sys.exit('Error: The grid file specified is either missing or lacking needed dimensions/variables.')

# ===================
# initial conditions
# ===================
# If starting from scratch, setup dimension variables and initial condition variables
if experiment in ('a', 'f'):
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

    # Assign initial condition variable values for EISMINT-2 experiment
    # Start with no ice
    thickness[:] = 0.0
    # zero velocity everywhere (only needed for HO solvers)
    normalVelocity[:] = 0.0
    # flat bed at sea level
    bedTopography[:] = 0.0
    # constant, arbitrary temperature, degrees C (doesn't matter since there is no ice initially)
    temperature[:] = 0.0 
    # Setup layerThicknessFractions
    layerThicknessFractions[:] = 1.0 / nVertLevels
# Now update/set origin location and distance array
r = ((xCell[:] - xsummit)**2 + (yCell[:] - ysummit)**2)**0.5

# ===================
# boundary conditions
# ===================
# Define values prescribed by Payne et al. 2000 paper.

params = exp_params[experiment]

# SMB field specified by EISMINT, constant in time for EISMINT-2
# It is a function of geographical position (not elevation)
Mmax = params['Mmax']  # [m/yr] maximum accumulation rate
Sb = params['Sb']      # gradient of accumulation rate change with horizontal distance  [m/a/km] 
Rel = params['Rel']    # [km]  accumulation rate 0 position
smb=numpy.minimum(Mmax, Sb * (Rel - r/1000.0)) # [m ice/yr]
SMB[:] = smb * rhoi / scyr  # in kg/m2/s

# Basal heat flux should be -4.2e-2
#G[:] = -4.2e-2  # [W/m2]

# Surface temperature
Tmin = params['Tmin'] # minimum surface air temperature [K]
ST = params['ST']     # gradient of air temperature change with horizontal distance [K/km]
#Tsfc[:] = Tmin + ST * r/1000.0

gridfile.close()
print 'Successfully added initial conditions for EISMINT2, experiment '+experiment+' to the file: ', filename

