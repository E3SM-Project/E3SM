#!/usr/bin/env python
'''
Plots profiles for hydro-ship test case C
This is the daily varying forcing.
It plots the daily time series for quantities of interest
for the last full some (?) days of the specified file.
It also plots the difference in the time series between the final two days.
'''
import sys
import numpy as np
import netCDF4
#import datetime
# import math
# from pylab import *
from optparse import OptionParser
import matplotlib.pyplot as plt
#from matplotlib import cm
# import time
#import random


parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file name to use", metavar="FILE")
#parser.add_option("-y", "--year", dest="year", type="int", help="year to visualize", metavar="YEAR")  # COULD MAKE A 'DAY' ARGUMENT
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")

options, args = parser.parse_args()

if not options.filename:
        options.filename = "output2.nc"
        print "No filename provided. Using:"+options.filename


class DailyData:

    def __init__(self, foutput, numDays):
        '''
        foutput = a netCDF4 object of the file to query
        numDays = number of days from the last day to assess
        '''
        self.numDays = numDays
        nCells = len(foutput.dimensions['nCells'])
        nEdges = len(foutput.dimensions['nEdges'])
        nTime = len(foutput.dimensions['Time'])
        xCell = foutput.variables['xCell'][:]
        yCell = foutput.variables['yCell'][:]
        yVertex = foutput.variables['yVertex'][:]
        xEdge = foutput.variables['xEdge'][:]
        #yEdge = f.variables['yEdge'][:]
        #areaCell = f.variables['areaCell'][:]
        #xtime= f.variables['xtime']
        days = foutput.variables['daysSinceStart'][:]
        days = days - days[0]  # relative to start time
        nt = 24  # hrs in a day
        self.hours = np.arange(24)

        h = f.variables['waterThickness']
        #u = f.variables['waterVelocityCellX']
        N = f.variables['effectivePressure']
        Q = f.variables['channelDischarge']

        indIce = np.where(h[0,:]>0.0)

        # identify needed days/indices
        lastFullDay = np.floor(days).max() - 1

        # calculate mean,min,max for all x values for needed variables
        self.N_mean = np.zeros((nt, numDays))
        self.N_min = np.zeros((nt, numDays))
        self.N_max = np.zeros((nt, numDays))
        self.h_mean = np.zeros((nt, numDays))
        self.h_min = np.zeros((nt, numDays))
        self.h_max = np.zeros((nt, numDays))
        self.Q_mean = np.zeros((nt, numDays))
        self.Q_max = np.zeros((nt, numDays))
        self.dayNum = np.zeros((numDays,))

        for i in range(numDays):
            thisDay = lastFullDay - i  # day number
            self.dayNum[-1-i] = thisDay  # Fill this from left to right.  Put last day on the left
            todayInd = np.where(np.logical_and(days >= thisDay, days < thisDay+1))[0]

            #print todayInd

            if len(todayInd) != 24:
               print "Error: Day {} does not have 24 time levels in it.".format(thisDay)
               print len(todayInd)

            for hr in range(len(todayInd)):  # loop through hours on this day
               N_now = N[todayInd[0]+hr, :]
               h_now = h[todayInd[0]+hr, :]
               Q_now = Q[todayInd[0]+hr, :]
               bigQ = np.where(Q_now>0.01)

               self.h_mean[hr,i] = h_now[indIce].mean()
               self.h_min[hr,i] = h_now[indIce].min()
               self.h_max[hr,i] = h_now[indIce].max()

               self.N_mean[hr,i] = N_now[indIce].mean()
               self.N_min[hr,i] = N_now[indIce].min()
               self.N_max[hr,i] = N_now[indIce].max()

               self.Q_max[hr,i] = Q_now.max()
               self.Q_mean[hr,i] = Q_now[bigQ].mean()


# Open file
f = netCDF4.Dataset(options.filename, 'r')
numDays = 5  # number of days before final too process
dayData = DailyData(f, numDays)  # analyze output data


# ===============
# set up figure 1 with empty axes
# left side axes show all days together
# ===============
fig = plt.figure(1, facecolor='w')
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(323, sharex=ax1)
ax3 = fig.add_subplot(325, sharex=ax1)


# Perform same plot operations for each day
for d in range(numDays):
   # water thickness mean/max
   plt.sca(ax1)
   #plt.plot(dayData.hours, dayData.h_mean[:,d], '-', label='h_mean{}'.format(dayData.dayNum[d]))
   #plt.plot(yr.days, yr.h_min, '--', color=yr.color, label='h_min'+yr.yearStr)
   plt.plot(dayData.hours, dayData.h_max[:,d], '--', label='h_max{}'.format(dayData.dayNum[d]))
   plt.ylabel('sheet thickness (m)')
   plt.title('Magnitude for all days')
   plt.legend()

   # N mean/max
   plt.sca(ax2)
   plt.plot(dayData.hours, dayData.N_mean[:,d]/1.0e6, '-', label='N_mean{}'.format(dayData.dayNum[d]))
   #plt.plot(days, N_min/1.0e6, 'b--', label='N_min')
   #plt.plot(days, N_max/1.0e6, 'b--', label='N_max')
   plt.ylabel('N (MPa)')
   plt.legend()

   # channel Q mean/max
   plt.sca(ax3)
   #plt.plot(dayData.hours, dayData.Q_mean[:,d], '-', label='Q_mean{}'.format(dayData.dayNum[d]))
   #plt.plot(days, Q_min, 'b--', label='h_min')
   plt.plot(dayData.hours, dayData.Q_max[:,d], '--', label='Q_max{}'.format(dayData.dayNum[d]))
   plt.xlabel('DOY')
   plt.ylabel('channel Q (m)')
   plt.legend()

plt.draw()
plt.show()
#  SCRIPT ONLY WORKING TO HERE

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

