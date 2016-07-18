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



f = netCDF4.Dataset(options.filename,'r')
#xtime = f.variables['xtime'][:]
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
#xEdge = f.variables['xEdge'][:]
#yEdge = f.variables['yEdge'][:]
h = f.variables['waterThickness'][time_slice,:]
#u = f.variables['waterVelocityCellX'][time_slice,:]
P = f.variables['waterPressure'][time_slice,:]
#H = f.variables['thickness'][time_slice,:]
opening = f.variables['openingRate'][time_slice,:]
closing = f.variables['closingRate'][time_slice,:]
melt = f.variables['basalMeltInput'][time_slice,:]
days = f.variables['daysSinceStart'][:]

print "Total number of time levels=", len(days)
print "Using time slice", time_slice, " which is year ", days[time_slice]/365.0

print "Attempting to read thickness field from landice_grid.nc."
fin = netCDF4.Dataset("landice_grid.nc",'r')
H = fin.variables['thickness'][0,:]


# Find center row  - currently files are set up to have central row at y=0
unique_ys=np.unique(yCell[:])
centerY=unique_ys[len(unique_ys)/2]
print "number of ys, center y index, center Y value", len(unique_ys), len(unique_ys)/2, centerY
ind = np.nonzero(yCell[:] == centerY)
x = xCell[ind]/1000.0

print "start plotting."

fig = plt.figure(1, facecolor='w')
# water thickness
ax = fig.add_subplot(211)
#plt.plot(x, H[ind]*917.0*9.81/1.0e5, '.-')
plt.plot(x, h[ind], '.-')
plt.xlabel('X-position (km)')
plt.ylabel('water depth (m)')
plt.grid(True)

# water pressure
ax = fig.add_subplot(212)
plt.plot(x, H[ind]*917.0*9.81, '.-')
plt.plot(x, P[ind], '.--')
plt.xlabel('X-position (km)')
plt.ylabel('water pressure (Pa)')
plt.grid(True)


# plot how close to SS we are
fig = plt.figure(2, facecolor='w')
ax1 = fig.add_subplot(311)
for i in ind:
    plt.plot(days/365.0, f.variables['waterThickness'][:,i])
plt.xlabel('Years since start')
plt.ylabel('water thickness (m)')
plt.grid(True)

ax = fig.add_subplot(312, sharex=ax1)
for i in ind:
    plt.plot(days/365.0, f.variables['effectivePressure'][:,i]/1.0e6)
plt.xlabel('Years since start')
plt.ylabel('effective pressure (MPa)')
plt.grid(True)

# plot opening/closing rates
ax = fig.add_subplot(313)
plt.plot(x, opening[ind], 'r', label='opening')
plt.plot(x, closing[ind], 'b', label='closing')
plt.plot(x, melt[ind] / 1000.0, 'g', label='melt')
plt.xlabel('X-position (km)')
plt.ylabel('rate (m/s)')
plt.legend()
plt.grid(True)




print "plotting complete"
plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

