#!/usr/bin/env python
# A script to compare MPAS model output to the EISMINT2 test cases.
# Matt Hoffman, LANL, December 2014

import sys, os
import datetime
try:
    import netCDF4
except ImportError:
    print 'Unable to import netCDF4 python modules:'
    sys.exit
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-e", "--experiment", dest="exp", type='string', help="experiment to setup", metavar="EXP")
#parser.add_option("-s", "--save", action="store_true", dest="saveimage", help="include this flag to save plot")
#parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plot")
options, args = parser.parse_args()

if not options.exp:
    sys.exit('Error: No experiment specified.  Please specify an experiment to visualize with the -e option')
experiment = options.exp.lower()  # allow lower or upper case to be input, but use lower case in the script logic
if experiment in ('a','b','c','d','f'):
    print 'Visualizing EISMINT2 Experiment ' + experiment
else:
    sys.exit("Error: Invalid experiment specified.  Please specify an experiment between 'a' and 'f', excluding 'd'")




################### DEFINE FUNCTIONS ######################
def xtime2numtime(xtime):
  """Define a function to convert xtime character array to numeric time values using datetime objects"""
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

def xtimeGetYear(xtime):
  """Get an array of years from an xtime array, ignoring any partial year information"""
  # First parse the xtime character array into a string
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function
  years = np.zeros( (len(xtimestr),) )
  for i in range(len(xtimestr)):
      years[i] = ( int(xtimestr[i].split('-')[0]) ) # Get the year part and make it an integer
  return years



def contourMPAS(field, contour_levs=None):
  """Contours irregular MPAS data on cells"""
  #-- Now let's grid your data.
  # First we'll make a regular grid to interpolate onto.
  numcols = nCells**0.5 * 4.0  # may want to adjust the density of the regular grid
  numrows = numcols
  xc = np.linspace(xCell.min(), xCell.max(), numcols)
  yc = np.linspace(yCell.min(), yCell.max(), numrows)
  xi, yi = np.meshgrid(xc, yc)
  #-- Interpolate at the points in xi, yi
  zi = griddata(xCell, yCell, field, xi, yi)
  #-- Display the results
  if contour_levs == None:
     im = plt.contour(xi, yi, zi)
  else:
     im = plt.contour(xi, yi, zi, contour_levs)

  #plt.scatter(xCell, yCell, c=temperature[timelev,:,-1], s=100, vmin=zi.min(), vmax=zi.max())  # to see the raw data on top
  plt.colorbar(im)
################### END OF FUNCTIONS ######################


# Find data for requested experiment
#experimentpath = "experiment_" + experiment + "/"
#filename = experimentpath + 'eismint2' + experiment + '.output.nc'
filename = "output.nc"

# open supplied MPAS output file and get variables needed
filein = netCDF4.Dataset(filename,'r')
xCell = filein.variables['xCell'][:]/1000.0
yCell = filein.variables['yCell'][:]/1000.0
xtime = filein.variables['xtime'][:]
nCells = len(filein.dimensions['nCells'])
nVertLevels = len(filein.dimensions['nVertLevels'])
#numtime = xtime2numtime(xtime)
years = xtimeGetYear(xtime)

thickness = filein.variables['thickness']
temperature = filein.variables['temperature']
basalTemperature = filein.variables['basalTemperature']
basalPmpTemperature = filein.variables['basalPmpTemperature']
flwa = filein.variables['flowParamA']
uReconstructX = filein.variables['uReconstructX']
uReconstructY = filein.variables['uReconstructY']
areaCell = filein.variables['areaCell'][:]
layerThicknessFractions = filein.variables['layerThicknessFractions'][:]

secInYr = 365.0*24.0*3600.0

timelev = -1  # Use final time
print 'Using final model time of ' + xtime[timelev,:].tostring().strip() + '\n'


##### Print some stats about the error
####print '===================================='
####print 'Max modeled thickness (m) = ' + str( thk.max() )
####print 'EISMINT models ice thickness at divide (m):'
####print '  3d models (10 of them): 2978.0 +/- 19.3'
####print '  2d models (3 of them):  2982.2 +/- 26.4'
####print '===================================='
####print ''

# ================
# ================
# Plot the results
# ================
# ================

# ================
# BASAL TEMPERATURE MAP
# ================

# make an educated guess about how big the markers should be.
if nCells**0.5 < 100.0:
  markersize=  max( int(round(  3600.0/(nCells**0.5)  )), 1)
  markershape = 'h'  # use hexes if the points are big enough, otherwise just dots
else:
  markersize=  max( int(round(  1800.0/(nCells**0.5)  )), 1)
  markershape = '.'
print 'Using a markersize of ', markersize


fig = plt.figure(1, facecolor='w')
fig.suptitle('Payne et al. Fig. 1, 3, 6, 9, or 11', fontsize=12, fontweight='bold')
#markersize = 30.0
#markershape='h'

# print ice locations with gray hexagons
iceIndices = np.where(thickness[timelev,:]>0.0)[0]
plt.scatter(xCell[iceIndices], yCell[iceIndices], markersize, (0.8, 0.8, 0.8), marker=markershape, edgecolors='none') # print ice locations with gray hexagons

# add contours of ice temperature over the top
contourMPAS(basalTemperature[timelev,:], np.linspace(240.0, 275.0, 8))

plt.axis('equal')
plt.title('Modeled basal temperature (K) \n at time ' + netCDF4.chartostring(xtime)[timelev].strip() )
plt.xlim( (0.0, 1500.0) ); plt.ylim( (0.0, 1500.0) )
plt.xlabel('X position (km)')
plt.ylabel('Y position (km)')




# ================
# STEADY STATE MAPS -  panels b and c are switched and with incorrect units in the paper
# ================
fig = plt.figure(2, facecolor='w', figsize=(12, 6), dpi=72)
fig.suptitle('Payne et al. Fig. 2 or 4', fontsize=12, fontweight='bold')

# ================
# panel a - thickness
ax1 = fig.add_subplot(131)

# print ice locations with gray hexagons
plt.scatter(xCell[iceIndices], yCell[iceIndices], markersize, (0.8, 0.8, 0.8), marker=markershape, edgecolors='none') # print ice locations with gray hexagons

# add contours of ice thickness over the top
contour_intervals = np.linspace(0.0, 5000.0,  5000.0/250.0+1)
contourMPAS(thickness[timelev,:], contour_levs=contour_intervals)
#contourMPAS(thickness[timelev,:])

plt.title('Final thickness (m)' )
plt.axis('equal')
#plt.xlim( (0.0, 750.0) ); plt.ylim( (0.0, 750.0) )
plt.xlabel('X position (km)'); plt.ylabel('Y position (km)')


# ================
# panel c - flux
ax = fig.add_subplot(133, sharex=ax1, sharey=ax1)

flux = np.zeros( (nCells,) )
for k in range(nVertLevels):
    speedLevel = (uReconstructX[timelev,:,k:k+2].mean(axis=1)**2 + uReconstructY[timelev,:,k:k+2].mean(axis=1)**2)**0.5
    flux += speedLevel * thickness[timelev,:] * layerThicknessFractions[k]

# print ice locations with gray hexagons
plt.scatter(xCell[iceIndices], yCell[iceIndices], markersize, (0.8, 0.8, 0.8), marker=markershape, edgecolors='none') # print ice locations with gray hexagons

# add contours over the top
#contour_intervals = np.linspace(0.0, 10.0,  10.0/0.5+1)
contour_intervals = np.linspace(0.0, 20.0,  11)
contourMPAS(flux * 3600.0*24.0*365.0 / 10000.0, contour_levs=contour_intervals)
#contourMPAS(flux * 3600.0*24.0*365.0 / 10000.0)
plt.axis('equal')
plt.title('Final flux (m$^2$ a$^{-1}$ / 10000)' )
#plt.xlim( (0.0, 750.0) ); plt.ylim( (0.0, 750.0) )
plt.xlabel('X position (km)'); plt.ylabel('Y position (km)')

# ================
# panel b - flow factor
ax = fig.add_subplot(132, sharex=ax1, sharey=ax1)

# print ice locations with gray hexagons
plt.scatter(xCell[iceIndices], yCell[iceIndices], markersize, (0.8, 0.8, 0.8), marker=markershape, edgecolors='none') # print ice locations with gray hexagons

# add contours over the top
#contour_intervals = np.linspace(0.0, 20.0, 20.0/2.0+1)
contour_intervals = np.linspace(0.0, 16.0, 16.0/0.5+1)
#contourMPAS(flwa[timelev,:,:].mean(axis=1) *3600.0*24.0*365.0 / 1.0e-17, contour_levs=contour_intervals)  # NOT SURE WHICH LEVEL FLWA SHOULD COME FROM - so taking column average
contourMPAS(flwa[timelev,:,:].mean(axis=1) * 3600.0*24.0*365.0 / 1.0e-17)  # NOT SURE WHICH LEVEL FLWA SHOULD COME FROM - so taking column average
plt.axis('equal')
plt.title('Final flow factor (10$^{-17}$ Pa$^{-3}$ a$^{-1}$)' )  # Note: the paper's figure claims units of 10$^{-25}$ Pa$^{-3}$ a$^{-1}$ but the time unit appears to be 10^-17
#plt.xlim( (0.0, 750.0) ); plt.ylim( (0.0, 750.0) )
plt.xlabel('X position (km)'); plt.ylabel('Y position (km)')



# ================
# DIVIDE EVOLUTION TIME SERIES
# ================
fig = plt.figure(3, facecolor='w')
fig.suptitle('Payne et al. Fig. 5, 7, or 8', fontsize=14, fontweight='bold')

# get indices for given time
if experiment =='b':
  endTime = 40000.0
else:
  endTime = 80000.0

# get index at divide - we set this up to be 750,750
divideIndex = np.logical_and( xCell == 750.0, yCell == 750.0)

# panel a - thickness
ax = fig.add_subplot(211)
timeInd = np.nonzero(years <= endTime)[0][0:] 
plt.plot( years[timeInd]/1000.0, thickness[timeInd, divideIndex], 'k.-')
plt.ylabel('Thickness (m)')

# panel b - basal temperature
ax = fig.add_subplot(212)
timeInd = np.nonzero(years <= endTime)[0][1:]  # skip the first index cause basalTemperature isn't calculated then
plt.plot( years[timeInd]/1000.0, basalTemperature[timeInd, divideIndex], 'k.-')
plt.ylabel('Basal temperature (K)')
plt.xlabel('Time (kyr)')




# ================
# TABLES
# ================
# Setup dictionaries of benchmark results for each experiment - values are mean, min, max from Tables in Payne et al. 2000
a_bench = {'stattype':'absolute', 'volume':(2.128, 2.060, 2.205), 'area':(1.034, 1.011, 1.097), 'meltfraction':(0.718, 0.587, 0.877), 'dividethickness':(3688.342, 3644.0, 3740.74), 'dividebasaltemp':(255.605, 254.16, 257.089)}
b_bench = {'stattype':'relative', 'volume':(-2.589, -3.079, -2.132), 'area':(0.0, 0.0, 0.0), 'meltfraction':(11.836, 3.307, 21.976), 'dividethickness':(-4.927, -5.387, -4.071), 'dividebasaltemp':(4.623, 4.47, 4.988)}
c_bench = {'stattype':'relative', 'volume':(-28.505, -29.226, -28.022), 'area':(-19.515, -20.369, -16.815), 'meltfraction':(-27.806, -39.353, -7.982), 'dividethickness':(-12.928, -13.948, -12.447), 'dividebasaltemp':(3.707, 3.389, 4.004)}
d_bench = {'stattype':'relative', 'volume':(-12.085, -12.890, -11.654), 'area':(-9.489, -10.184, -6.924), 'meltfraction':(-1.613, -4.744, 1.001), 'dividethickness':(-2.181, -2.517, -1.985), 'dividebasaltemp':(-0.188, -0.209, -0.149)}
f_bench = {'stattype':'absolute', 'volume':(0.0, 0.0, 0.0), 'area':(0.0, 0.0, 0.0), 'meltfraction':(0.0, 0.0, 0.0), 'dividethickness':(0.0, 0.0, 0.0), 'dividebasaltemp':(0.0, 0.0, 0.0)}
# Setup dictionary of dictionaries for each experiment
benchmarks = {'a':a_bench, 'b':b_bench, 'c':c_bench, 'd':d_bench, 'f':f_bench}


bench = benchmarks[experiment]  # Get the benchmark dictionary

fig = plt.figure(4, facecolor='w')
fig.suptitle('Payne et al. Table 4, 5, 6, 7, 8, or 9: showing min/mean/max of community', fontsize=12, fontweight='bold')

fig.add_subplot(151)
volume = (thickness[timelev,iceIndices] * areaCell[iceIndices]).sum() / 1000.0**3 / 10.0**6
plt.plot( np.zeros((3,)), bench['volume'], 'k*')  # benchmark results
if bench['stattype'] == 'relative':
    initIceIndices = np.where(thickness[0,:]>0.0)[0]
    volume = (volume / ((thickness[0,initIceIndices] * areaCell[initIceIndices]).sum() / 1000.0**3 / 10.0**6) - 1.0) * 100.0
    plt.ylabel('Volume change (%)')
else:
    plt.ylabel('Volume (10$^6$ km$^3$)')
plt.plot( (0.0,), volume, 'ro')  # MPAS results
plt.xticks(())

fig.add_subplot(152)
area = (areaCell[iceIndices]).sum() / 1000.0**2 / 10.0**6
areaAbsolute = area
plt.plot( np.zeros((3,)), bench['area'], 'k*')  # benchmark results
if bench['stattype'] == 'relative':
    initArea = (areaCell[initIceIndices]).sum() / 1000.0**2 / 10.0**6
    area = (area / initArea - 1.0) * 100.0
    plt.ylabel('Area change (%)')
else:
    plt.ylabel('Area (10$^6$ km$^2$)')
plt.plot( (0.0,), area, 'ro')  # MPAS results
plt.xticks(())

fig.add_subplot(153)
warmBedIndices = np.where(np.logical_and(thickness[timelev,:] > 0.0, basalTemperature[timelev,:] >= (basalPmpTemperature[timelev,:] - 0.01) ) )[0]  # using threshold here to identify melted locations
meltfraction = areaCell[warmBedIndices].sum() / 1000.0**2 / 10.0**6 / areaAbsolute
plt.plot( np.zeros((3,)), bench['meltfraction'], 'k*')  # benchmark results
if bench['stattype'] == 'relative':
    # use time 1 instead of 0 since these fields aren't fully populated at time 0
    initIceIndices = np.where(thickness[1,:]>0.0)[0]
    initArea = (areaCell[initIceIndices].sum() / 1000.0**2 / 10.0**6)
    initWarmBedIndices = np.where(np.logical_and(thickness[1,:] > 0.0, basalTemperature[1,:] >= (basalPmpTemperature[1,:] - 0.01) ) )[0]  # using threshold here to identify melted locations
    initWarmArea = (areaCell[initWarmBedIndices].sum() / 1000.0**2 / 10.0**6)
    initMeltFraction = initWarmArea / initArea
    meltfraction = (meltfraction / initMeltFraction - 1.0) * 100.0
    plt.ylabel('Melt fraction change (%)')
else:
    plt.ylabel('Melt fraction')
plt.plot( (0.0,), meltfraction, 'ro')  # MPAS results
plt.xticks(())

fig.add_subplot(154)
dividethickness = thickness[timelev, divideIndex]
plt.plot( np.zeros((3,)), bench['dividethickness'], 'k*')  # benchmark results
if bench['stattype'] == 'relative':
    dividethickness = (dividethickness / thickness[0, divideIndex] - 1.0) * 100.0
    plt.ylabel('Divide thickness change (%)')
else:
    plt.ylabel('Divide thickness (m)')
plt.plot( (0.0,), dividethickness, 'ro')  # MPAS results
plt.xticks(())

fig.add_subplot(155)
dividebasaltemp = basalTemperature[timelev, divideIndex]
plt.plot( np.zeros((3,)), bench['dividebasaltemp'], 'k*')  # benchmark results
if bench['stattype'] == 'relative':
    # use time 1 instead of 0 since these fields aren't fully populated at time 0
    dividebasaltemp = dividebasaltemp - basalTemperature[1, divideIndex]
    plt.ylabel('Divide basal temp. change (K)')
else:
    plt.ylabel('Divide basal temp. (K)')
plt.plot( (0.0,), dividebasaltemp, 'ro')  # MPAS results
plt.xticks(())

plt.tight_layout()





plt.draw()
plt.show()


#if options.saveimage:
#    plotname = experimentpath + 'EISMINT2-'+experiment+'-basaltemp.png'
#    plt.savefig(plotname, dpi=150)
#    print 'Saved plot as ' + plotname

#if options.hidefigs:
#     print "Plot display disabled with -n argument."
#else:
#     print 'Showing plot...  Close plot window to exit.'
#     plt.show()


