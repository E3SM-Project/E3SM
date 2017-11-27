#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
'''

import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm

outfname = 'globalStats.nc'
runs=[ adir for adir in sorted(os.listdir('.')) if (os.path.isdir(adir) and os.path.isfile(os.path.join(adir, outfname)))]
print "Original run list:", runs

# reorder to put the 'control' runs at the beginning
special_runs = ('steady', 'no-melt')
for r in special_runs:
   if r in runs:
      runs.remove(r)
      runs.insert(0, r)

# optionally exclude some subset
#runs[:] = [r for r in runs if not 'amp300_per02' in r]

print "Will process the following directories: ", runs

# --- Define some needed parameters
rhoi = 910.0  # ice density
windowLength = 7      # window length for smoothing some time series
# ----------------

# ------ Needed functions ------

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


# --- Define data structures ---

#class variabilitySet:

class modelRun:
   def __init__(self, run):
      '''
      This reads results from a model run and saves and analyzes the needed results.

      run = name of subdir in which the run was performed
      '''
      f = netCDF4.Dataset(run + '/' + outfname, 'r')
      self.nt = len(f.dimensions['Time'])
      self.yrs = np.zeros((self.nt,))
      #yrs = f.variables['daysSinceStart'][:] / 365.0
      self.xtimes = f.variables['xtime'][:]
   
   
      self.yrs = xtime2numtimeMy(self.xtimes)
      self.dyrs = self.yrs[1:] - self.yrs[0:-1]
   
      # resampled version of time array - needed for filtering.
      # (Filtering needed b/c the occasional tiny time step that the model is exhibiting leads to inaccurate (noisy) derivatives)
      self.resampYrs = np.linspace(0.0, self.yrs.max(), len(self.yrs)*2)
      self.dresampYrs = self.resampYrs[1:] - self.resampYrs[0:-1]
   
      self.melt = f.variables['avgSubshelfMelt'][:]
      self.totalmelt = f.variables['totalFloatingBasalMassBal'][:] / 1.0e12
      self.cumumelt = self.totalmelt*0.0
      for t in np.arange(self.nt-1)+1:
         self.cumumelt[t] = self.cumumelt[t-1] + self.totalmelt[t-1] * self.dyrs[t-1]
   
      self.VAF = f.variables['volumeAboveFloatation'][:] / 1.0e12 * rhoi
      self.resampVAF = np.interp(self.resampYrs, self.yrs, self.VAF) # generate y values for each x
      self.VAFsmooth = scipy.signal.savgol_filter(self.resampVAF, window_length=windowLength, polyorder=1)
      #VAFrate = (VAF[1:]-VAF[0:-1] ) / dyrs
      self.VAFsmoothrate = (self.VAFsmooth[1:] - self.VAFsmooth[0:-1] ) / self.dresampYrs
   
      self.GA = f.variables['groundedIceArea'][:] / 1.0e3**2
      self.resampGA = np.interp(self.resampYrs, self.yrs, self.GA) # generate y values for each x
      self.GAsmooth = scipy.signal.savgol_filter(self.resampGA, window_length=windowLength, polyorder=1)
      self.GArate = (self.GAsmooth[1:] - self.GAsmooth[0:-1]) / self.dresampYrs


      f.close()



# --- set up figure axes ---

print "Setting up figure axes."

fig = plt.figure(1, facecolor='w')

nrow=4
ncol=2

xtickSpacing = 20.0

# melt forcing
axMeanMelt = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('mean melt (m/yr)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axTotalMelt = fig.add_subplot(nrow, ncol, 3, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('total melt (Gt/yr)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axCumuMelt = fig.add_subplot(nrow, ncol, 5, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('cumulative melt (Gt)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

# VAF
axVAF = fig.add_subplot(nrow, ncol, 2, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF (Gt)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axVAFrate = fig.add_subplot(nrow, ncol, 4, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF rate (Gt/yr)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

# grounded area
axGA = fig.add_subplot(nrow, ncol, 6, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('Grounded area (km^2)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

axGArate = fig.add_subplot(nrow, ncol, 8, sharex=axMeanMelt)
plt.xlabel('Year')
plt.ylabel('Grounded area rate (km^2/yr)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()



# =================== repeat cleaner for presentation ===========

fig2 = plt.figure(2, facecolor='w')

nrow=3
ncol=1

# melt forcing
ax2MeanMelt = fig2.add_subplot(nrow, ncol, 1)
plt.xlabel('Year')
plt.ylabel('mean melt (m/yr)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

# VAF
ax2VAF = fig2.add_subplot(nrow, ncol, 2, sharex=ax2MeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF (Gt)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()

ax2VAFrate = fig2.add_subplot(nrow, ncol, 3, sharex=ax2MeanMelt)
plt.xlabel('Year')
plt.ylabel('VAF rate (Gt/yr)')
plt.xticks(np.arange(22)*xtickSpacing)
plt.grid()


# =========
print "Done setting up figure axes."
# =========

# --- Define colors for lines ---
#colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan']
n150 = sum("amp150" in r for r in runs)
colors = [ cm.autumn(x) for x in np.linspace(0.0, 1.0, n150) ]
n300 = sum("amp300" in r for r in runs)
colors.extend( [ cm.winter(x) for x in np.linspace(0.0, 1.0, n300) ] )
color_index = 0


# ================
# Loop over runs and collect needed data
# ================
runData = {}  # init empty dictionary
for run in runs:
   print "Processing run: " + run

   # Build a list that contains all run data objects
   runData[run] = modelRun(run)

print "Processing complete.\n"

# ================
# Loop over runs and plot data
# ================
for run in runs:
   print "Plotting run: " + run

   thisRun = runData[run]
   # Pull needed data for plotting from this run
   # TODO: replace local variables below in plotting commands with reference to object variables
   yrs=thisRun.yrs
   melt=thisRun.melt
   totalmelt=thisRun.totalmelt
   cumumelt=thisRun.cumumelt
   VAF=thisRun.VAF
   resampYrs=thisRun.resampYrs
   VAFsmoothrate=thisRun.VAFsmoothrate
   GA=thisRun.GA
   GArate=thisRun.GArate

   # optional rules about coloring
   if "steady" in run:
      lw = 3
      color = 'k'
      #colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, len(runs)) ];  color = colors[color_index];  color_index += 1; lw=1
   elif run=="no-melt":
      lw = 2
      color = [0.7, 0.7, 0.7]  # gray
   else:
      lw = 1
      color = colors[color_index]
      color_index += 1


   # actual plotting begins here
   axMeanMelt.plot(yrs, melt, color=color, linewidth=lw)
   axTotalMelt.plot(yrs, totalmelt, color=color, linewidth=lw)
   axCumuMelt.plot(yrs, cumumelt, color=color, linewidth=lw)

   axVAF.plot(yrs, VAF, label = run, color=color, linewidth=lw)
   #axVAFrate.plot(yrs[1:], VAFrate, color=color, linewidth=lw)
   axVAFrate.plot(resampYrs[1:], VAFsmoothrate, color=color, linewidth=lw)

   axGA.plot(yrs, GA, color=color, linewidth=lw)
   axGArate.plot(resampYrs[1:], GArate, color=color, linewidth=lw)


   # actual plotting begins here for second figure
   ax2MeanMelt.plot(yrs, melt, color=color, linewidth=lw)
   ax2VAF.plot(yrs, VAF, label = run, color=color, linewidth=lw)
   #ax2VAFrate.plot(yrs[1:], VAFrate, color=color, linewidth=lw)
   #ax2VAFrate.plot(resampYrs[1:], VAFsmoothrate, color=color, linewidth=lw)
   ax2VAFrate.plot(resampYrs[windowLength:], VAFsmoothrate[windowLength-1:], color=color, linewidth=lw)
   #ax2VAFrate.plot(resampYrs[:windowLength], VAFsmoothrate[0:windowLength-1], 'x', color=color)  # used to see what values we are throwing away



axVAF.legend(loc='best', ncol=2)

# second ticks for SLR
def GTtoSL(GT):
   #return (GT-VAF[0])*1.0e12/1028.0/1000.0**3/362.0e6*1000.0*1000.0
   return GT *1.0e12/1028.0/1000.0**3/362.0e6*1000.0*1000.0 * -1.0

axSLR=axVAF.twinx()
y1, y2=axVAF.get_ylim()
x1, x2=axVAF.get_xlim()
axSLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLR.set_ylabel('S.L. equiv. (mm)')
axSLR.set_xlim(x1, x2)

axSLRrate=axVAFrate.twinx()
y1, y2=axVAFrate.get_ylim()
x1, x2=axVAFrate.get_xlim()
axSLRrate.set_ylim(GTtoSL(y1), GTtoSL(y2))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
axSLRrate.set_ylabel('S.L. equiv. (mm/yr)')
axSLRrate.set_xlim(x1, x2)





# ===================  for figure 2
handles, labels = ax2VAF.get_legend_handles_labels()
l1 = ax2VAF.legend(handles[0:2], labels[0:2], loc='lower left', ncol=1)  # control runs
#l2 = ax2VAF.legend(handles[2:], labels[2:], loc='lower center', ncol=2)  # variability runs
ax2VAF.add_artist(l1)

ax2SLR=ax2VAF.twinx()
y1, y2=ax2VAF.get_ylim()
x1, x2=ax2VAF.get_xlim()
ax2SLR.set_ylim(GTtoSL(y1) - GTtoSL(VAF[0]), GTtoSL(y2) - GTtoSL(VAF[0]))
#axSLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
ax2SLR.set_ylabel('S.L. equiv. (mm)')
ax2SLR.set_xlim(x1, x2)

ax2SLRrate=ax2VAFrate.twinx()
y1, y2=ax2VAFrate.get_ylim()
x1, x2=ax2VAFrate.get_xlim()
ax2SLRrate.set_ylim(GTtoSL(y1), GTtoSL(y2))
#ax2SLR.set_yticks( range(int(GTtoSL(y1)), int(GTtoSL(y2))) )
ax2SLRrate.set_ylabel('S.L. equiv. (mm/yr)')
ax2SLRrate.set_xlim(x1, x2)



plt.show()
