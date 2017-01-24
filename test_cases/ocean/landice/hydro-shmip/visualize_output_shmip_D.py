#!/usr/bin/env python
'''
Plots profiles for hydro-ship test case D
This is the seasonally varying forcing.
It plots the annual time series for quantities of interest
for the year specified along with the previous year (assumed to be
in files each one year long).
It also plots the difference in the time series between the two years.
'''
import sys
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
import random

secInYr = 3600.0 * 24.0 * 365.0  # Note: this may be slightly wrong for some calendar types!

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="template for file to visualize of format: FILE.YYYY.nc", metavar="FILE")
parser.add_option("-y", "--year", dest="year", type="int", help="year to visualize", metavar="YEAR")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using 'output2'."
        options.filename = "output2"

if not options.year:
	print "No year provided. Using year 0001."
        options.year = 1

class AnnualData:

    def __init__(self, year):
        self.year = year
        self.yearStr = '{:04d}'.format(year)
        self.color = 'b'
        self.filename = options.filename + "." + self.yearStr + ".nc"
        f = netCDF4.Dataset(self.filename, 'r')   
        nCells = len(f.dimensions['nCells'])
        nEdges = len(f.dimensions['nEdges'])
        nTime = len(f.dimensions['Time'])
        xCell = f.variables['xCell'][:]
        yCell = f.variables['yCell'][:]
        yVertex = f.variables['yVertex'][:]
        xEdge = f.variables['xEdge'][:]
        #yEdge = f.variables['yEdge'][:]
        #areaCell = f.variables['areaCell'][:]
        #xtime= f.variables['xtime']
        days = f.variables['daysSinceStart'][:]
        self.days = days - days[0]  # relative to start time
        
        # ---------------
        # Plot 1: time series over the year of various quantities
        # ---------------
        
        
        h = f.variables['waterThickness']
        #u = f.variables['waterVelocityCellX']
        N = f.variables['effectivePressure']
        Q = f.variables['channelDischarge']
        
        indIce = np.where(h[0,:]>0.0)
        
        # calculate mean,min,max for all x values for needed variables
        self.N_mean = np.zeros((nTime,))
        self.N_min = np.zeros((nTime,))
        self.N_max = np.zeros((nTime,))
        self.h_mean = np.zeros((nTime,))
        self.h_min = np.zeros((nTime,))
        self.h_max = np.zeros((nTime,))
        self.Q_mean = np.zeros((nTime,))
        self.Q_max = np.zeros((nTime,))
        for i in range(nTime):
            N_now = N[i,:]
            h_now = h[i,:]
            Q_now = Q[i,:]
            bigQ = np.where(Q_now>0.01)
        
            self.h_mean[i] = h_now[indIce].mean()
            self.h_min[i] = h_now[indIce].min()
            self.h_max[i] = h_now[indIce].max()
        
            self.N_mean[i] = N_now[indIce].mean()
            self.N_min[i] = N_now[indIce].min()
            self.N_max[i] = N_now[indIce].max()
        
            self.Q_max[i] = Q_now.max()
            self.Q_mean[i] = Q_now[bigQ].mean()
      
     
# Create objects for each year
thisyear = options.year
thisYearData = AnnualData(thisyear)
thisYearData.color = 'b'

lastyear = thisyear-1
lastYearData = AnnualData(lastyear)
lastYearData.color = 'c'

# list of year objects to iterate over
years = (lastYearData, thisYearData)

# ===============
# set up figure 1 with empty axes
# left side axes show both years together
# ===============
fig = plt.figure(1, facecolor='w')
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(323, sharex=ax1)
ax3 = fig.add_subplot(325, sharex=ax1)

#years = (lastyear, thisyear)

# Perform same plot operations for each year
for yr in years:
   # water thickness mean/max
   plt.sca(ax1)
   plt.plot(yr.days, yr.h_mean, '-', color=yr.color, label='h_mean'+yr.yearStr)
   #plt.plot(yr.days, yr.h_min, '--', color=yr.color, label='h_min'+yr.yearStr)
   plt.plot(yr.days, yr.h_max, '--', color=yr.color, label='h_max'+yr.yearStr)
   plt.ylabel('sheet thickness (m)')
   plt.title('Magnitude for both years')
   plt.legend()
   
   # N mean/max
   plt.sca(ax2)
   plt.plot(yr.days, yr.N_mean/1.0e6, '-', color=yr.color, label='N_mean'+yr.yearStr)
   #plt.plot(days, N_min/1.0e6, 'b--', label='N_min')
   #plt.plot(days, N_max/1.0e6, 'b--', label='N_max')
   plt.ylabel('N (MPa)')
   plt.legend()
   
   # channel Q mean/max
   plt.sca(ax3)
   plt.plot(yr.days, yr.Q_mean, '-', color=yr.color, label='Q_mean'+yr.yearStr)
   #plt.plot(days, Q_min, 'b--', label='h_min')
   plt.plot(yr.days, yr.Q_max, '--', color=yr.color, label='Q_max'+yr.yearStr)
   plt.xlabel('DOY')
   plt.ylabel('channel Q (m)')
   plt.legend()

#plt.draw()


# ===============
# Plot differnces on axes on right side of figure
# ===============
if 1: #(thisYearData.days == lastYearData.days).all():
   ax4 = fig.add_subplot(322, sharex=ax1)
   plt.plot(thisYearData.days, thisYearData.h_mean - lastYearData.h_mean, 'r-', label='h_mean diff')
   plt.plot(thisYearData.days, thisYearData.h_max - lastYearData.h_max, 'r--', label='h_max diff')
   plt.title('Difference between years')
   plt.legend()

   ax5 = fig.add_subplot(324, sharex=ax1)
   plt.plot(thisYearData.days, (thisYearData.N_mean - lastYearData.N_mean)/1.0e6, '-', label='N_mean diff')
   plt.legend()

   ax6 = fig.add_subplot(326, sharex=ax1)
   plt.plot(thisYearData.days, thisYearData.Q_mean - lastYearData.Q_mean, 'r-', label='Q_mean diff')
   plt.plot(thisYearData.days, thisYearData.Q_max - lastYearData.Q_max, 'r--', label='Q_max diff')
   plt.legend()
   plt.xlabel('DOY')

else:
   print "Time indices for the two years do not match. Skipping difference plot."
   print "This year times:", thisYearData.days
   print "Last year times:", lastYearData.days



plt.draw()

print "plotting complete"

if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')

if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

