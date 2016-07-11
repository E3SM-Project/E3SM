#!/usr/bin/env python
'''
Plots profiles for hydro-margin test case
'''
import numpy as np
import netCDF4
#import datetime
# import math
# from pylab import *
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib import cm
# from matplotlib.contour import QuadContourSet
# import time

secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize", metavar="FILE")
parser.add_option("-t", "--time", dest="time", help="time step to visualize (0 based)", metavar="TIME")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

if not options.time:
	print "No time provided. Using time -1."
        time_slice = -1
else:
        time_slice = int(options.time)


fig = plt.figure(1, facecolor='w')

f = netCDF4.Dataset(options.filename,'r')

#xtime = f.variables['xtime'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
#xEdge = f.variables['xEdge'][:]
#yEdge = f.variables['yEdge'][:]
h = f.variables['waterThickness'][time_slice,:]
u = f.variables['waterVelocityCellX'][time_slice,:]
N = f.variables['effectivePressure'][time_slice,:]


# Find center row  - currently files are set up to have central row at y=0
unique_ys=np.unique(yCell[:])
centerY=unique_ys[len(unique_ys)/2]
print "number of ys, center y index, center Y value", len(unique_ys), len(unique_ys)/2, centerY
ind = np.nonzero(yCell[:] == centerY)
x = xCell[ind]/1000.0

print "start plotting."

ax = fig.add_subplot(311)
plt.plot(x, h[ind] * u[ind] * 10000.0, '.-')
plt.xlabel('X-position (km)')
plt.ylabel('water flux (m^3/s)')


ax = fig.add_subplot(312)
plt.plot(x, h[ind]*10000.0, '.-')
plt.xlabel('X-position (km)')
plt.ylabel('xsect area (m^2)')
plt.yscale('log')

ax = fig.add_subplot(313)
plt.plot(x, N[ind], '.-')
plt.xlabel('X-position (km)')
plt.ylabel('effective pressure (Pa)')

print "plotting complete"

plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

