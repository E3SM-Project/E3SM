#!/usr/bin/env python
'''
Script to compare some scalar values from different runs of Thwaites melt variability experiment.
'''


# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-p", dest="process", action="store_true", help="read and process data", metavar="FILE")
options, args = parser.parse_args()


import sys
import os
import netCDF4
import datetime
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from matplotlib import cm
import pickle


outfname = 'output.nc'
runs=[ adir for adir in sorted(os.listdir('.')) if (os.path.isdir(adir) and os.path.isfile(os.path.join(adir, outfname)))]
print "Original run list:", runs
#runs = [ 'steady', 'amp300_per20_pha0.00','amp300_per20_pha0.25', 'amp300_per20_pha0.50', 'amp300_per20_pha0.75']
runs = ['steady','amp300_per20_pha0.00']
runs = ['steady']
runs = ['amp300_per20_pha0.00']
# reorder to put the 'control' runs at the beginning
special_runs = ('steady', 'no-melt')
for r in special_runs:
   if r in runs:
      runs.remove(r)
      runs.insert(0, r)
print "Will process the following directories: ", runs

# ----- needed functions ----
def xtime2numtimeMy(xtime):
  """Define a function to convert xtime character array to numeric time values using local arithmetic"""
  # First parse the xtime character array into a string
  xtimestr = netCDF4.chartostring(xtime) # convert from the character array to an array of strings using the netCDF4 module's function

  numtime = np.zeros( (len(xtimestr),) )
  ii = 0
  dayOfMonthStart = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
  for stritem in xtimestr:
      itemarray = stritem.strip().replace('_', '-').replace(':', '-').split('-')  # Get an array of strings that are Y,M,D,h,m,s
      results = [int(i) for i in itemarray]
      numtime[ii] = results[0] + (dayOfMonthStart[results[1]-1]-1 + results[2]) / 365.0  # decimal year
      ii += 1
  return numtime


# --- Define data structures ---

class modelRun:
   def __init__(self, run):
      '''
      This reads results from a model run and saves and analyzes the needed results.

      run = name of subdir in which the run was performed
      '''

      # some metadata about run
      if 'amp' in run:
         self.amp = float(run[3:6])
         self.per = float(run[10:12])
         self.pha = float(run[16:])
      else:
         self.amp = 0.0
         self.per = 0.0
         self.pha = 0.0

      #print self.amp, self.per, self.pha

      f = netCDF4.Dataset(run + '/' + outfname, 'r')
      self.nt = len(f.dimensions['Time'])
      self.yrs = np.zeros((self.nt,))
      #yrs = f.variables['daysSinceStart'][:] / 365.0
      self.xtimes = f.variables['xtime'][:]


      self.yrs = xtime2numtimeMy(self.xtimes)
      self.dyrs = self.yrs[1:] - self.yrs[0:-1]


      xCell = f.variables['xCell'][:]
      yCell = f.variables['yCell'][:]
      xEdge  = f.variables['xEdge'][:]
      yEdge  = f.variables['yEdge'][:]
      conc= f.variables['cellsOnCell'][:,:]
      cone= f.variables['cellsOnEdge'][:,:]
      dcEdge = f.variables['dcEdge'][:]
      areaCell = f.variables['areaCell'][:]
      topg = f.variables['bedTopography'][0,:]

      nt = len(f.dimensions['Time'])

      self.time = np.zeros((nt,))
      self.mnGLslope = np.zeros((nt,))
      self.GA = np.zeros((nt,))
      self.mnGLmelt = np.zeros((nt,))
      self.mnGLdepthMPASMethod = np.zeros((nt,))
      self.mnGLdepth = np.zeros((nt,))
      self.mnCFmelt = np.zeros((nt,))
      self.mnGLbedslope = np.zeros((nt,))
      self.mnGLxpos = np.zeros((nt,))

      spdThresh = 100.0/3.14e7
      spdThresh = 900.0/3.14e7
      for t in range(nt):
#      for t in np.arange(0,nt,10):
#      for t in np.arange(0,284,1):
         print t
         usrf = f.variables['upperSurface'][t,:]
         edgeMask = f.variables['edgeMask'][t,:]
         cellMask = f.variables['cellMask'][t,:]
         thk = f.variables['thickness'][t,:]
         speed = f.variables['surfaceSpeed'][t,:]
         melt = f.variables['floatingBasalMassBal'][t,:]

         ind = np.where(edgeMask&256>0)[0]
         skipped = 0
         kept = 0
         for i in range(len(ind)):
             c1 = cone[ind[i]-1,0]-1 # in python 0-based 
             c2 = cone[ind[i]-1,1]-1 # in python 0-based
             if thk[c1] <= 0.0 or thk[c2] <= 0.0 or c1 <= 0 or c2 <= 0 or speed[c1]<spdThresh or speed[c2]<spdThresh:
                skipped += 1
                pass
             kept += 1
             slp = np.absolute( (usrf[c1] - usrf[c2]) / dcEdge[ind[i]-1] )
             self.mnGLslope[t] += slp # accumulate for mean calc below
             self.mnGLmelt[t] += -1.0 * min(melt[c1], melt[c2]) # get the melt rate on the floating side of the GL (the grounded side will be 0)
             self.mnGLdepth[t] += (topg[c1] + topg[c2])/2.0 # approx GL depth here by avg of two neighboring cells TODO: improve this

             self.mnGLxpos[t] += xEdge[i]
             # bed slope at GL
             if (cellMask[c1] & 256 > 0):
                sign = 1.0
             else:
                sign = -1.0
             self.mnGLbedslope[t] += sign*(topg[c1] - topg[c2])/dcEdge[ind[i]-1]
         #print "skipped {} of {} edges".format(skipped, len(ind))
         self.mnGLslope[t] /= kept  # calc mean
         self.mnGLmelt[t] /= kept
         self.mnGLdepth[t] /= kept
         self.mnGLbedslope[t] /= kept
         self.mnGLxpos[t] /= kept
         self.GA[t] = (areaCell[:] * (910.0/1028.0 * thk > -1.0 * topg)).sum()

         ind = np.where(np.logical_and(edgeMask&256>0, np.logical_and(yEdge<-400.0e3, yEdge>-500.0e3)))[0]
         self.mnGLxpos[t] = xEdge[ind].mean()
         # calc GL depth using MPAS method
         ind = np.where(cellMask&256>0)[0]
         frac = 0.25
         frac = 0.5
         GLdepths = np.sort(topg[ind])
         #print GLdepths
         fracIdx = int(round(len(GLdepths) * frac))  # get the index that represents the deepest x% of entries
         self.mnGLdepthMPASMethod[t] = GLdepths[0:fracIdx].mean()

         # CF melt
         ind = np.where(xCell<-1.56e6)
         self.mnCFmelt[t] = melt[ind].mean()

# ================
# Loop over runs and collect needed data
# ================
runData = {}  # init empty dictionary
groups = {} # groups are a dictionary of dictionaries.  May want to changes this to a dictionary of custom objects that can include metadata about the group alongside the group's data children.  KISS for now.
for run in runs:
   print "Processing run: " + run

   # Build a list that contains all run data objects
   if options.process:
      runData[run] = modelRun(run)
      file_pi = open(run + '/' +'calc_gl_data.pi', 'w') 
      pickle.dump(runData[run], file_pi)
      file_pi.close()
   else:
      file_pi = open(run + '/' +'calc_gl_data.pi', 'r') 
      runData[run] = pickle.load(file_pi)
      file_pi.close()

   # Populate each run into its group for ensemble runs (not for controls)
   if 'amp' in run:
      groupName = run[0:12]  # this is the amp and per
      #print groupName
      phase = run[13:]
      #print phase
      if groupName not in groups:
         groups[groupName] = {} # init an empty dict to add group members
      groups[groupName][phase] = runData[run]  # stick the run object into this group

print "Processing complete.\n"


# ------- figure --------
fig = plt.figure(1, facecolor='w')

ax1 = fig.add_subplot(1, 3, 1)
plt.xlabel('time')
plt.ylabel('GA')
plt.grid()


ax2 = fig.add_subplot(1, 3, 2)
plt.xlabel('time')
plt.ylabel('slope')
plt.grid()


axGAslope = fig.add_subplot(1, 3, 3)
plt.xlabel('GA')
plt.ylabel('slope')
plt.grid()

# ------
fig = plt.figure(2, facecolor='w')
axGAslope2 = fig.add_subplot(1, 1, 1)
plt.xlabel('change in GA from initial')
plt.ylabel('mean slope at GL')
plt.grid()


# -----
fig = plt.figure(3, facecolor='w')
nrow = 2; ncol = 2
axGLslope = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('time (yr)')
plt.ylabel('mean slope at GL')
plt.grid()

axGLmelt = fig.add_subplot(nrow, ncol, 2)
plt.xlabel('time (yr)')
plt.ylabel('mean melt at GL')
plt.grid()

axGLdepth = fig.add_subplot(nrow, ncol, 3)
plt.xlabel('time (yr)')
plt.ylabel('mean depth at GL')
plt.grid()

axGAvsGLmelt = fig.add_subplot(nrow, ncol, 4)
plt.xlabel('change in GA from initial')
plt.ylabel('melt at GL')
plt.grid()


fig = plt.figure(4, facecolor='w')
nrow = 2; ncol = 2
axCFmelt = fig.add_subplot(nrow, ncol, 1)
plt.xlabel('time (yr)')
plt.ylabel('mean melt of CF')
plt.grid()

axBedslope = fig.add_subplot(nrow, ncol, 2)
plt.xlabel('time (yr)')
plt.ylabel('mean bed slope at GL')
plt.grid()

axGLrate= fig.add_subplot(nrow, ncol, 3)
plt.xlabel('time (yr)')
plt.ylabel('GL xpos retreat rate (m/yr)')
plt.grid()


# --- Define colors for lines ---
colors = []
#colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan']
n150 = sum("amp150" in r for r in runs)
colors.extend( [ cm.autumn(x) for x in np.linspace(0.0, 1.0, n150) ])
n300 = sum("amp300" in r for r in runs)
colors.extend( [ cm.winter(x) for x in np.linspace(0.0, 1.0, n300) ] )
color_index = 0


# ================
# Loop over runs and plot data
# ================
runNumber = 0
for run in runs:
   print "Plotting run: " + run

   thisRun = runData[run]
   stuff =  thisRun.mnGLdepthMPASMethod # somehow i need to access this member before using...?
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


   ax1.plot(thisRun.GA, '.', color=color)
   ax2.plot(thisRun.mnGLslope, '.', color=color)
   axGAslope.plot(thisRun.GA, thisRun.mnGLslope, '.', color=color)
   axGAslope2.plot(thisRun.GA[0]-thisRun.GA, thisRun.mnGLslope, '.', color=color, markersize=2)
   #fig 3
   axGLslope.plot(thisRun.yrs, thisRun.mnGLslope, color=color)
   axGLmelt.plot(thisRun.yrs, thisRun.mnGLmelt, color=color)
   axGLdepth.plot(thisRun.yrs, thisRun.mnGLdepth, color=color)
   axGLdepth.plot(thisRun.yrs, thisRun.mnGLdepthMPASMethod[:], '.', color=color)
   step = 40
   axGLdepth.plot(thisRun.yrs[step/2:-1*step/2], 20.0*(thisRun.mnGLdepthMPASMethod[:-1*step]-thisRun.mnGLdepthMPASMethod[step:]) / (thisRun.yrs[:-1*step]-thisRun.yrs[step:]), color='g', lineWidth=1.5)
   step = 50
   axGLdepth.plot(thisRun.yrs[step/2:-1*step/2], 80.0*(thisRun.mnGLdepth[:-1*step]-thisRun.mnGLdepth[step:]) / (thisRun.yrs[:-1*step]-thisRun.yrs[step:]), color='m', lineWidth=1.5)

   axGAvsGLmelt.plot(thisRun.GA[0]-thisRun.GA, thisRun.mnGLmelt, '.', color=color, markersize=1)
   axCFmelt.plot(thisRun.yrs, thisRun.mnCFmelt, color=color)
   axBedslope.plot(thisRun.yrs, thisRun.mnGLbedslope, color=color)
   #axGLrate.plot(thisRun.yrs[1:], -1*(thisRun.mnGLxpos[1:]-thisRun.mnGLxpos[:-1]), color=color, lineWidth=0.5)
   step=30
   axGLrate.plot(thisRun.yrs[step/2:-1*step/2], -1.0*(thisRun.mnGLxpos[:-1*step]-thisRun.mnGLxpos[step:]) / (thisRun.yrs[:-1*step]-thisRun.yrs[step:]), color=color, lineWidth=1.5)


plt.show()
