#!/usr/bin/env python
# A script to compare MPAS model output to the Halfar analytic solution of the dome test case.
# Matt Hoffman, LANL, September 2013

import sys
import datetime
try:
    import netCDF4
except ImportError:
    print 'Unable to import netCDF4 python modules:'
    sys.exit
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to test", metavar="FILE")
parser.add_option("-t", "--time", dest="t", help="which time level to use", metavar="T")
parser.add_option("-s", "--save", action="store_true", dest="saveimage", help="include this flag to save plot")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plot")

options, args = parser.parse_args()
if not options.filename:
   options.filename = 'output.nc'
   print 'No file specified.  Attempting to use output.nc'
if options.t:
   timelev = int(options.t)
else:
   timelev = -1
   print 'No time level specified.  Attempting to use final time.'

################### DEFINE FUNCTIONS ######################
# Define the function to calculate the Halfar thickness
def halfar(t,x,y, A, n, rho):
  # A   # s^{-1} Pa^{-3}
  # n   # Glen flow law exponent
  # rho # ice density kg m^{-3}

  # These constants should come from setup_dome_initial_conditions.py.
  # For now they are being hardcoded.
  R0 = 60000.0 * np.sqrt(0.125)   # initial dome radius
  H0 = 2000.0 * np.sqrt(0.125)    # initial dome thickness at center
  g = 9.80616  # gravity m/s/s -- value used by mpas_constants
  alpha = 1.0/9.0
  beta = 1.0/18.0
  secpera = 3600.0*24.0*365.0  # value for gregorian_noleap calendar
  Gamma = 2.0 / (n+2.0) * A * (rho * g)**n

  x0 = 0.0; y0 = 0.0   # The IC file puts the center of the dome and the domain at (0,0)

  t0 = (beta/Gamma) * (7.0/4.0)**3 * (R0**4/H0**7)  # NOTE: These constants assume n=3 - they need to be generalized to allow other n's 
  t=t+t0
  t=t/t0

  H=np.zeros(len(x))
  for i in range(len(x)):
      r = np.sqrt( (x[i] - x0)**2 + (y[i] - y0)**2)
      r=r/R0
      inside = max(0.0, 1.0 - (r / t**beta)**((n+1.0) / n))

      H[i] = H0 * inside**(n / (2.0*n+1.0)) / t**alpha
  return H

# Define a function to convert xtime character array to numeric time values using datetime objects
def xtime2numtime(xtime):
  # First parse the xtime character array into a string 
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function

  dt = []
  for stritem in xtimestr:
      itemarray = stritem.strip().replace('_', '-').replace(':', '-').split('-')  # Get an array of strings that are Y,M,D,h,m,s
      results = [int(i) for i in itemarray]
      if (results[0] < 1900):  # datetime has a bug where years less than 1900 are invalid on some systems
         results[0] += 1900
      dt.append( datetime.datetime(*results) ) # * notation passes in the array as arguments

  numtime = netCDF4.date2num(dt, units='seconds since '+str(dt[0]))   # use the netCDF4 module's function for converting a datetime to a time number
  return numtime

################### END OF FUNCTIONS ######################


# open supplied MPAS output file and get thickness slice needed
filein = netCDF4.Dataset(options.filename,'r')
xCell = filein.variables['xCell'][:]
yCell = filein.variables['yCell'][:]
xtime = filein.variables['xtime'][:]

thk = filein.variables['thickness'][:]
xtime = filein.variables['xtime'][:] 
numtime = xtime2numtime(xtime)

# Find out what the ice density and flowA values for this run were.
print '\nCollecting parameter values from the output file.'

flowA = filein.config_default_flowParamA
print 'Using a flowParamA value of: ' + str(flowA)
flow_n = filein.config_flowLawExponent
print 'Using a flowLawExponent value of: ' + str(flow_n)
if flow_n != 3:
        print 'Error: The Halfar script currently only supports a flow law exponent of 3.'
        sys.exit
rhoi = filein.config_ice_density
print 'Using an ice density value of: ' + str(rhoi)

dynamicThickness = filein.config_dynamic_thickness
print '\nDynamic thickness for this run = ' + str(dynamicThickness)

print 'Using model time of ' + xtime[timelev,:].tostring().strip() + '\n'

if filein.config_calendar_type != "gregorian_noleap":
        print 'Error: The Halfar script currently assumes a gregorian_noleap calendar.  Modify it to proceed with your calendar type of: ', filein.config_calendar_type
        sys.exit


# Call the halfar function
thkHalfar = halfar(numtime[timelev]-numtime[0], xCell, yCell, flowA, flow_n, rhoi)

thkDiff = (thk[timelev, :] - thkHalfar) 
thkDiffIce = thkDiff[ np.where( thk[timelev,:] > 0.0) ]  # Restrict to cells modeled to have ice
RMS = ( (thkDiffIce**2).sum() / float(len(thkDiffIce)) )**0.5


# Print some stats about the error
print 'Error statistics for cells modeled to have ice:'
print '* RMS error = ' + str( RMS )
print '* Minimum error = ' + str( thkDiffIce.min() )
print '* Maximum error = ' + str( thkDiffIce.max() )
print '* Mean error = ' + str( thkDiffIce.mean() )
print '* Median error = ' + str( np.median(thkDiffIce) )
print '* Mean absolute error = ' + str( np.absolute(thkDiffIce).mean() )
print '* Median absolute error = ' + str( np.median(np.absolute(thkDiffIce)) )
print ''

# Plot the results
fig = plt.figure(1, figsize=(16, 4.5), facecolor='w', dpi=100)
markersize = 35.0
gray = np.ones(3)*0.8

fig.add_subplot(1,3,1)
maskindices = np.nonzero(thk[:][timelev,:] > 0.0)[:]
plt.scatter(xCell/1000.0,yCell/1000.0,markersize,gray, marker='.', edgecolors='none')
plt.scatter(xCell[maskindices]/1000.0,yCell[maskindices]/1000.0,markersize,thk[timelev,maskindices], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness (m) \n at time ' + netCDF4.chartostring(xtime)[timelev].strip() ) 
plt.xlabel('x (km)'); plt.ylabel('y (km)')

fig.add_subplot(1,3,2)
halmaskindices = np.nonzero(thkHalfar[:] > 0.0)[:]
plt.scatter(xCell/1000.0,yCell/1000.0,markersize,gray, marker='.', edgecolors='none')
plt.scatter(xCell[halmaskindices]/1000.0,yCell[halmaskindices]/1000.0,markersize,thkHalfar[halmaskindices], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Analytic thickness (m) \n at time ' + netCDF4.chartostring(xtime)[timelev].strip() ) 
plt.xlabel('x (km)'); plt.ylabel('y (km)')

fig.add_subplot(1,3,3)
plt.scatter(xCell/1000.0,yCell/1000.0,markersize,gray, marker='.', edgecolors='none')
plt.scatter(xCell/1000.0,yCell/1000.0,markersize,thkDiff, marker='h', edgecolors='none', cmap=plt.cm.RdBu)
plt.colorbar()
plt.clim([-1.0*np.absolute(thkDiff).max(), np.absolute(thkDiff).max()])
plt.axis('equal')
plt.title('Modeled thickness - Analytic thickness (m) \n at time ' + netCDF4.chartostring(xtime)[timelev].strip() ) 
plt.xlabel('x (km)'); plt.ylabel('y (km)')

plt.draw()

if options.saveimage:
    plotname = 'halfar-results.png'
    plt.savefig(plotname, dpi=150)
    print 'Saved plot as ' + plotname

if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     print 'Showing plot...  Close plot window to exit.'
     plt.show()


