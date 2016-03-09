#!/usr/bin/env python
'''
Plots velocity profiles for a diagnostic solve for a range of resolutions, with and without GLP.
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

reslist = (10000, 5000, 2000, 1000, 500, 250)
GLbit = 256
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
	print "No time provided. Using time 0."
        time_slice = 0
else:
        time_slice = int(options.time)


################### DEFINE FUNCTIONS ######################
def get_data(filename):

  f = netCDF4.Dataset(filename,'r')

  #xtime = f.variables['xtime'][:]
  xCell = f.variables['xCell'][:]
  yCell = f.variables['yCell'][:]
  #xEdge = f.variables['xEdge'][:]
  #yEdge = f.variables['yEdge'][:]
  surfaceSpeed = f.variables['surfaceSpeed'][time_slice,:]
  #edgeMask = f.variables['edgeMask']  # just get the object

  #vert_levs = len(f.dimensions['nVertLevels'])

  # Find center row  - currently files are set up to have central row at y=0
  ind = np.nonzero( yCell == 0.0 )
  x = xCell[ind]/1000.0
  u = surfaceSpeed[ind]*secInYr

  f.close()
  return x, u

colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, len(reslist)) ]
fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111)
for i in range(len(reslist)):
  res = reslist[i]
  # no glp first
  fname = "{}m.nc".format(res)
  print "Processing file", fname
  x, u = get_data(fname)
  plt.plot(x, u, '.-', color=colors[i], label="{}m, no GLP".format(res))
  # glp next
  fname = "{}m-glp.nc".format(res)
  print "Processing file", fname
  x, u = get_data(fname)
  plt.plot(x, u, '.--', color=colors[i], label="{}m, GLP".format(res))


plt.xlabel('X-position (km)')
plt.ylabel('Speed (m/yr)')
plt.title('Profile at y=0')
plt.legend()
plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:     
     plt.show()

