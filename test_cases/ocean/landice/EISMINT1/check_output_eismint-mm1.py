#!/usr/bin/env python
# A script to compare MPAS model output to the EISMINT Moving Margin 1 test case.
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
#numtime = xtime2numtime(xtime)

# Find out what the ice density and flowA values for this run were.
print '\nCollecting parameter values from the output file.'

flowA = filein.config_default_flowParamA
print 'Using a flowParamA value of: ' + str(flowA)
flow_n = filein.config_flowLawExponent
print 'Using a flowLawExponent value of: ' + str(flow_n)
rhoi = filein.config_ice_density
print 'Using an ice density value of: ' + str(rhoi)
dynamicThickness = filein.config_dynamic_thickness
print 'Dynamic thickness for this run = ' + str(dynamicThickness)

print 'Using model time of ' + xtime[timelev,:].tostring().strip() + '\n'


# Print some stats about the error
print '===================================='
print 'Max modeled thickness (m) = ' + str( thk.max() )
print 'EISMINT models ice thickness at divide (m):'
print '  3d models (10 of them): 2978.0 +/- 19.3'
print '  2d models (3 of them):  2982.2 +/- 26.4'
print '===================================='
print ''

# Plot the results
fig = plt.figure(1, facecolor='w')
markersize = 30.0

fig.add_subplot(1,1,1)
plt.scatter(xCell,yCell,markersize,thk[timelev,:], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness (m) \n at time ' + netCDF4.chartostring(xtime)[timelev].strip() ) 


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


