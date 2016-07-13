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
yVertex = f.variables['yVertex'][:]
#xEdge = f.variables['xEdge'][:]
#yEdge = f.variables['yEdge'][:]
h = f.variables['waterThickness'][time_slice,:]
u = f.variables['waterVelocityCellX'][time_slice,:]
N = f.variables['effectivePressure'][time_slice,:]
days = f.variables['daysSinceStart'][:]
basalmelt = f.variables['basalMeltInput'][time_slice,:]
surfmelt = f.variables['externalWaterInput'][time_slice,:]
xtime= f.variables['xtime']
areaCell = f.variables['areaCell'][:]


#print "attempting to get input data from landice_grid.nc!"
#fin = netCDF4.Dataset('landice_grid.nc','r')
#H = fin.variables['thickness'][0,:]

print "Using time level ", time_slice, ", which is xtime=",''.join( xtime[time_slice,:])


totalMelt = ( (basalmelt + surfmelt)/1000.0 * areaCell * (h>0.0)).sum() / (yVertex.max() - yVertex.min()) / 2.0 * 10000.0  # 10000.0 is the 10 km width that Hewitt uses; 2 is because we want half the domain since we reflect it
print "Total water input volume rate for half the domain (m^3/s)=", totalMelt

# Find center row  - currently files are set up to have central row at y=0
unique_ys=np.unique(yCell[:])
centerY=unique_ys[len(unique_ys)/2]
print "number of ys, center y index, center Y value", len(unique_ys), len(unique_ys)/2, centerY
ind = np.nonzero(yCell[:] == centerY)
x = xCell[ind]/1000.0

print "start plotting."

fig = plt.figure(1, facecolor='w')
ax1 = fig.add_subplot(311)
plt.plot(x, h[ind] * u[ind] * 10000.0, '.-')
plt.plot(1000.0, totalMelt, 'r*')
plt.xlabel('X-position (km)')
plt.ylabel('water flux (m^3/s)')
plt.grid(True)


ax = fig.add_subplot(312, sharex=ax1)
plt.plot(x, h[ind]*10000.0, '.-')
plt.xlabel('X-position (km)')
plt.ylabel('xsect area (m^2)')
plt.yscale('log')
plt.grid(True)

ax = fig.add_subplot(313, sharex=ax1)
plt.plot(x, N[ind]/1.0e6, '.-', label='modeled transient to SS')
plt.plot(x, (basalmelt[ind]/900.0*1.0e13/h[ind]) / 1.0e6, '.--r', label='SS N=f(h)')  # steady state N=f(h) from the cavity evolution eqn
plt.xlabel('X-position (km)')
plt.ylabel('effective pressure (MPa)')
plt.ylim((0.0, 0.1))
plt.grid(True)
plt.legend()



# plot how close to SS we are
fig = plt.figure(2, facecolor='w')
ax1 = fig.add_subplot(211)
for i in ind:
    plt.plot(days/365.0, f.variables['waterThickness'][:,i])
plt.xlabel('Years since start')
plt.ylabel('water thickness (m)')
plt.grid(True)

ax = fig.add_subplot(212, sharex=ax1)
for i in ind:
    plt.plot(days/365.0, f.variables['effectivePressure'][:,i]/1.0e6)
plt.xlabel('Years since start')
plt.ylabel('effective pressure (MPa)')
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

