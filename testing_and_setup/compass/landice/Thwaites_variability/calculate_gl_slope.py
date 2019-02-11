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

# reorder to put the 'control' runs at the beginning
special_runs = ('steady', 'no-melt')
for r in special_runs:
   if r in runs:
      runs.remove(r)
#      runs.insert(0, r)
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
      conc= f.variables['cellsOnCell'][:,:]
      cone= f.variables['cellsOnEdge'][:,:]
      dcEdge = f.variables['dcEdge'][:]
      areaCell = f.variables['areaCell'][:]
      topg = f.variables['bedTopography'][0,:]

      nt = len(f.dimensions['Time'])

      self.mnGLslope = np.zeros((nt,))
      self.GA = np.zeros((nt,))
      self.mnGLmelt = np.zeros((nt,))
      self.mnGLdepthMPASMethod = np.zeros((nt,))
      self.mnGLdepth = np.zeros((nt,))
      self.time = np.zeros((nt,))
      for t in range(nt):
      #for t in np.arange(0,nt,10):
         print t
         usrf = f.variables['upperSurface'][t,:]
         edgeMask = f.variables['edgeMask'][t,:]
         cellMask = f.variables['cellMask'][t,:]
         thk = f.variables['thickness'][t,:]
         melt = f.variables['floatingBasalMassBal'][t,:]
         ind = np.where(edgeMask&256>0)[0]
         for i in range(len(ind)):
             c1 = cone[ind[i]-1,0]-1 # in python 0-based 
             c2 = cone[ind[i]-1,1]-1 # in python 0-based
             slp = np.absolute( (usrf[c1] - usrf[c2]) / dcEdge[ind[i]-1] )
             self.mnGLslope[t] += slp # accumulate for mean calc below
             self.mnGLmelt[t] += -1.0 * min(melt[c1], melt[c2]) # get the melt rate on the floating side of the GL (the grounded side will be 0)
             self.mnGLdepth[t] += (topg[c1] + topg[c2])/2.0 # approx GL depth here by avg of two neighboring cells TODO: improve this
         self.mnGLslope[t] /= len(ind)  # calc mean
         self.mnGLmelt[t] /= len(ind)
         self.mnGLdepth[t] /= len(ind)
         self.GA[t] = (areaCell[:] * (910.0/1028.0 * thk > -1.0 * topg)).sum()

         # calc GL depth using MPAS method
         ind = np.where(cellMask&256>0)[0]
         frac = 0.25
         GLdepths = np.sort(topg[ind])
         #print GLdepths
         fracIdx = int(round(len(GLdepths) * frac))  # get the index that represents the deepest x% of entries
         self.mnGLdepthMPASMethod[t] = GLdepths[0:fracIdx].mean()

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


# --- Define colors for lines ---
colors = []
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan']
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

plt.show()
