#!/usr/bin/env python
import numpy as np
import netCDF4
import datetime
# import math
# from pylab import *
from optparse import OptionParser
import matplotlib.pyplot as plt
import fnmatch, os
# from matplotlib.contour import QuadContourSet
# import time

template = "output-mask*.nc"
GLbit = 256
secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to analyze.  If not provided, the script will look for a filename template of the form: '{}'".format(template), metavar="FILE")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")
# parser.add_option("-v", "--var", dest="variable", help="variable to visualize", metavar="VAR")
# parser.add_option("--max", dest="maximum", help="maximum for color bar", metavar="MAX")
# parser.add_option("--min", dest="minimum", help="minimum for color bar", metavar="MIN")

options, args = parser.parse_args()

if options.filename:
        print "Using single file:", options.filename
        filelist = (options.filename,)
else:
	print "No filename provided. Using template."
        filelist = fnmatch.filter(os.listdir('.'), template)
        filelist.sort()  # sort this cause the list will be in arbitrary order
print "List of files to process:", filelist

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
  numtime /= (3600.0 * 24.0 * 365.0)
  numtime -= numtime[0]  # return years from start
  return numtime

def xtimeGetYear(xtime):
  """Get an array of years from an xtime array, ignoring any partial year information"""
  # First parse the xtime character array into a string 
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function
  years = np.zeros( (len(xtimestr),) )
  for i in range(len(xtimestr)):
      years[i] = ( int(xtimestr[i].split('-')[0]) ) # Get the year part and make it an integer
  return years



# These are only needed if gridding.
# nx = 30
# ny = 35

GLposAll = np.zeros((1,4))

for filename in filelist:
  print "Processing file:", filename

  f = netCDF4.Dataset(filename,'r')
  
  xtime = f.variables['xtime'][:]
  #years = xtime2numtime(xtime)
  years = xtimeGetYear(xtime)
  #thickness = f.variables['thickness'][time_slice,:]
  #dcedge = f.variables['dcEdge'][:]
  #bedTopography = f.variables['bedTopography']  # not needed
  #xCell = f.variables['xCell'][:]
  #yCell = f.variables['yCell'][:]
  xEdge = f.variables['xEdge'][:]
  yEdge = f.variables['yEdge'][:]
  #angleEdge = f.variables['angleEdge'][:]
  #lowerSurface = f.variables['lowerSurface'][time_slice,:]
  #upperSurface = f.variables['upperSurface'][time_slice,:]
  #surfaceSpeed = f.variables['surfaceSpeed'][time_slice,:]
  #basalSpeed = f.variables['basalSpeed'][time_slice,:]
  #floatingEdges = f.variables['floatingEdges'][time_slice,:]
  edgeMask = f.variables['edgeMask']  # just get the object
  #normalVelocity = f.variables['normalVelocity']
  #uReconstructX = f.variables['uReconstructX']
  #uReconstructX = f.variables['uReconstructX']
  #uReconstructY = f.variables['uReconstructY']
  
  vert_levs = len(f.dimensions['nVertLevels'])
  nt = len(f.dimensions['Time'])
  
  # print "nx = ", nx, " ny = ", ny
  print "vert_levs = ", vert_levs, " time_length = ", nt
  
  GLpos = np.zeros((nt,4))  # time, min, mean, max.  This array is just this file.
  for t in range(nt):
    GLpos[t,0] = years[t]
    GLind = np.nonzero( np.logical_and( ( (edgeMask[t,:] & GLbit) / GLbit == 1), (xEdge > 0.0) ) )
    #print 'Time, GL position values', years[t], xEdge[GLind]
    if len(GLind[0]) > 0:
      GLpos[t,1] = xEdge[GLind].min() / 1000.0
      GLpos[t,2] = xEdge[GLind].mean() / 1000.0
      GLpos[t,3] = xEdge[GLind].max() / 1000.0
    else:
      GLpos[t,1:] = 0.0
  GLposAll = np.append(GLposAll, GLpos, axis=0)  # add on this file
  f.close()


fig = plt.figure(1, facecolor='w')
ax = fig.add_subplot(111)
# Calculate GL position
print "Final GL position (time, min, mean, max):", GLpos[-1,:]

plt.plot(GLposAll[:,0], GLposAll[:,1], ':b')
plt.plot(GLposAll[:,0], GLposAll[:,2], '-bo')
plt.plot(GLposAll[:,0], GLposAll[:,3], ':b')
plt.xlabel('Time (yrs from start)')
plt.ylabel('GL position (km)')
plt.title('GL position over time. final pos={:.1f}'.format(GLpos[-1,2]))


ax2 = ax.twinx()
ax2.plot(GLposAll[:,0], GLposAll[:,3] - GLposAll[:,1], '-r')
ax2.set_ylabel('GL range (km)', color='r')
for tl in ax2.get_yticklabels():
    tl.set_color('r')


plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:     
     plt.show()


