#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
'''

import sys
import os
import netCDF4
import datetime
from math import sqrt
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import scipy.signal

#runs=[
#'no-melt',
#'amp0',
#'amp200per2',
#'amp200per20',
#'amp200per70',
#'amp500per2',
#'amp500per20',
#'amp500per70'
#]
outfname = 'globalStats.nc.clean'
runs=[ adir for adir in sorted(os.listdir('.')) if (os.path.isdir(adir) and os.path.isfile(os.path.join(adir, outfname)))]

print "Will process the following directories: ", runs

fig = plt.figure(1, facecolor='w')

nrow=4
ncol=2

# melt forcing
axMeanMelt = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('mean melt (m/yr)')
plt.xticks(np.arange(22)*20.0)
plt.grid()

axTotalMelt = fig.add_subplot(nrow, ncol, 3, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('total melt (Gt/yr)')
plt.xticks(np.arange(22)*20.0)
plt.grid()

axCumuMelt = fig.add_subplot(nrow, ncol, 5, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('cumulative melt (Gt)')
plt.xticks(np.arange(22)*20.0)
plt.grid()

# VAF
axVAF = fig.add_subplot(nrow, ncol, 2, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF (Gt)')
plt.xticks(np.arange(22)*20.0)
plt.grid()

axVAFrate = fig.add_subplot(nrow, ncol, 4, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF rate (Gt/yr)')
plt.xticks(np.arange(22)*20.0)
plt.grid()

# grounded area
axGA = fig.add_subplot(nrow, ncol, 6, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('Grounded area (km^2)')
plt.xticks(np.arange(22)*20.0)
plt.grid()

axGArate = fig.add_subplot(nrow, ncol, 8, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('Grounded area rate (km^2/yr)')
plt.xticks(np.arange(22)*20.0)
plt.grid()



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
  return numtime / (3600.0*24.0*365.0) # in years

def xtime2numtimeMy(xtime):
  """Define a function to convert xtime character array to numeric time values using local arithmetic"""
  # First parse the xtime character array into a string
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function

  numtime = np.zeros( (len(xtimestr),) )
  ii = 0
  for stritem in xtimestr:
      itemarray = stritem.strip().replace('_', '-').replace(':', '-').split('-')  # Get an array of strings that are Y,M,D,h,m,s
      results = [int(i) for i in itemarray]
      numtime[ii] = results[0] + results[1]/12.0 + results[2]/365.0  # approximate years
      ii += 1
  return numtime


rhoi = 910.0

for run in runs:
   print "Processing run: " + run
   f = netCDF4.Dataset(run + '/' + outfname, 'r')
   nt = len(f.dimensions['Time'])
   yrs = np.zeros((nt,))
   #yrs = f.variables['daysSinceStart'][:] / 365.0
   xtimes = f.variables['xtime'][:]

   yrs = xtime2numtimeMy(xtimes)
   dyrs = yrs[1:]-yrs[0:-1]
   melt = f.variables['avgSubshelfMelt'][:]
   totalmelt = f.variables['totalFloatingBasalMassBal'][:] / 1.0e12
   cumumelt = totalmelt*0.0
   for t in np.arange(nt-1)+1:
      cumumelt[t] = cumumelt[t-1] + totalmelt[t-1] * dyrs[t-1]

   VAF = f.variables['volumeAboveFloatation'][:] / 1.0e12 * rhoi
   VAFrate = (VAF[1:]-VAF[0:-1] ) / dyrs

   GA = f.variables['groundedIceArea'][:] / 1.0e3**2
   GAsmooth = scipy.signal.savgol_filter(GA, window_length=27, polyorder=2)
   GArate = (GAsmooth[1:] - GAsmooth[0:-1]) / dyrs

   axMeanMelt.plot(yrs, melt)
   axTotalMelt.plot(yrs, totalmelt)
   axCumuMelt.plot(yrs, cumumelt)

   axVAF.plot(yrs, VAF, label = run)
   axVAFrate.plot(yrs[1:], VAFrate)

   axGA.plot(yrs, GA)
   axGArate.plot(yrs[1:], GArate)

   f.close()

axVAF.legend(loc='best')

# second ticks for SLR
def GTtoSL(GT):
   return (GT-VAF[0])*1.0e12/1028.0/1000.0**3/362.0e6*1000.0*1000.0

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1), GTtoSL(y2))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('sea level equivalent (mm)')
axSLR.set_xlim(x1, x2)

plt.show()
