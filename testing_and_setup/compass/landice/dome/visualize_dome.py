#!/usr/bin/env python
import numpy
# from netCDF import *
# import math
from Scientific.IO.NetCDF import *
# from pylab import *
from optparse import OptionParser
import matplotlib.pyplot as plt
# from matplotlib.contour import QuadContourSet
# import time


parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize", metavar="FILE")
parser.add_option("-t", "--time", dest="time", help="time step to visualize (0 based)", metavar="TIME")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")
# parser.add_option("-v", "--var", dest="variable", help="variable to visualize", metavar="VAR")
# parser.add_option("--max", dest="maximum", help="maximum for color bar", metavar="MAX")
# parser.add_option("--min", dest="minimum", help="minimum for color bar", metavar="MIN")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "output.nc"

if not options.time:
	print "No time provided. Using time 0."
        time_slice = 0
else:
        time_slice = int(options.time)

#if not options.variable:
#	parser.error("Variable is a required input.")

# if not options.maximum:
#      	color_max = 0.0
# else:
# 	color_max = float(options.maximum)

# if not options.minimum:
# 	color_min = 0.0
# else:
# 	color_min = float(options.minimum)

# These are only needed if gridding.
# nx = 30
# ny = 35


secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!


f = NetCDFFile(options.filename,'r')

times = f.variables['xtime']
thickness = f.variables['thickness']
dcedge = f.variables['dcEdge']
#bedTopography = f.variables['bedTopography']  # not needed
xCell = f.variables['xCell']
yCell = f.variables['yCell']
xEdge = f.variables['xEdge']
yEdge = f.variables['yEdge']
angleEdge = f.variables['angleEdge']
temperature = f.variables['temperature']
lowerSurface = f.variables['lowerSurface']
upperSurface = f.variables['upperSurface']
normalVelocity = f.variables['normalVelocity']
uReconstructX = f.variables['uReconstructX']
uReconstructX = f.variables['uReconstructX']
uReconstructY = f.variables['uReconstructY']

vert_levs = f.dimensions['nVertLevels']

time_length = times.shape[0]

# print "nx = ", nx, " ny = ", ny
print "vert_levs = ", vert_levs, " time_length = ", time_length


# print "Computing global max and min"
# junk = thickness[:,:,0]
# maxval = junk.max()
# minval = junk.min()
# 
# del junk
# 
# junk = thickness[:,:]
# global_max = junk.max()
# global_min = junk.min()


var_slice = thickness[time_slice,:]
# var_slice = var_slice.reshape(time_length, ny, nx)


# print "Global max = ", global_max, " Global min = ", global_min
# print "Surface max = ", maxval, " Surface min = ", minval

fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111, aspect='equal')
# C = plt.contourf(xCell, yCell, var_slice )
plt.scatter(xCell[:], yCell[:], 80, var_slice, marker='h', edgecolors='none')
plt.colorbar()
plt.title('thickness at time ' + str(time_slice) )
plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('dome_thickness.png')

fig = plt.figure(2)
ax = fig.add_subplot(121, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, lowerSurface[time_slice,:], marker='h', edgecolors='none')
plt.colorbar()
plt.title('lower surface at time ' + str(time_slice) )
plt.draw()
ax = fig.add_subplot(122, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, upperSurface[time_slice,:], marker='h', edgecolors='none')
plt.colorbar()
plt.title('upper surface at time ' + str(time_slice) )
plt.draw()
if options.saveimages:
        plt.savefig('dome_surfaces.png')


fig = plt.figure(3)
for templevel in xrange(0,vert_levs):
    ax = fig.add_subplot(3,4,templevel+1, aspect='equal')
    var_slice = temperature[time_slice,:,templevel]
    # C = plt.contourf(xCell, yCell, var_slice )
    plt.scatter(xCell[:], yCell[:], 40, var_slice, marker='h', edgecolors='none')
    plt.colorbar()
    plt.title('temperature at level '+ str(templevel) + ' at time ' + str(time_slice) )
    plt.draw()
if options.saveimages:
        plt.savefig('dome_temperature.png')

fig = plt.figure(4)
ax = fig.add_subplot(121, aspect='equal')
plt.scatter(xEdge[:], yEdge[:], 80, normalVelocity[time_slice,:,vert_levs-1] * secInYr, marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xEdge[:], yEdge[:], numpy.cos(angleEdge[:]) * normalVelocity[time_slice,:,vert_levs-1] * secInYr, numpy.sin(angleEdge[:]) * normalVelocity[time_slice,:,vert_levs-1] * secInYr)
plt.title('normalVelocity of bottom layer at time ' + str(time_slice) )
plt.draw()
ax = fig.add_subplot(122, aspect='equal')
plt.scatter(xEdge[:], yEdge[:], 80, normalVelocity[time_slice,:,0] * secInYr, marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xEdge[:], yEdge[:], numpy.cos(angleEdge[:]) * normalVelocity[time_slice,:,0] * secInYr, numpy.sin(angleEdge[:]) * normalVelocity[time_slice,:, 0] * secInYr )
plt.title('normalVelocity of top layer at time ' + str(time_slice) )
plt.draw()
if options.saveimages:
        plt.savefig('dome_normalVelocity.png')

fig = plt.figure(5, facecolor='w')
ax = fig.add_subplot(121, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, uReconstructX[time_slice,:,0] * secInYr, marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xCell[:], yCell[:], uReconstructX[time_slice,:, 0] * secInYr, uReconstructY[time_slice,:, 0] * secInYr )
plt.title('uReconstructX of top layer at time ' + str(time_slice) )
plt.draw()
ax = fig.add_subplot(122, aspect='equal')
plt.scatter(xCell[:], yCell[:], 80, uReconstructY[time_slice,:,0] * secInYr, marker='h', edgecolors='none')
plt.colorbar()
plt.quiver(xCell[:], yCell[:], uReconstructX[time_slice,:, 0] * secInYr, uReconstructY[time_slice,:, 0] * secInYr )
plt.title('uReconstructY of top layer at time ' + str(time_slice) )
plt.draw()
if options.saveimages:
        plt.savefig('dome_uReconstruct.png')



if options.hidefigs:
     print "Plot display disabled with -n argument."
else:     
     plt.show()

f.close()

