#!/usr/bin/env python
import numpy
from netCDF4 import Dataset

from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import time

import glob


parser = OptionParser()
           
options, args = parser.parse_args()

iteration = int(args[0])

varNames = ['kineticEnergyCellMax', 'kineticEnergyCellAvg','layerThicknessMin','temperatureAvg']
timeUnit = 'hours'
daysPerUnit = 1/24.

inFolder = '.'
inFiles = glob.glob('%s/analysis_members/global*.nc'%(inFolder))
if(len(inFiles) == 0):
  print "Error: files not found in %s"%(inFolder)
  exit(1)
inFiles.sort()

fields = []
for varIndex in range(len(varNames)):
  fields.append(numpy.empty(0))
times = numpy.empty(0)
for inFileName in inFiles:
  inFile = Dataset(inFileName,'r')

  localTime = inFile.variables['daysSinceStartOfSim']
  for varIndex in range(len(varNames)):
    fieldLocal = numpy.array(inFile.variables[varNames[varIndex]])
    fields[varIndex] = numpy.append(fields[varIndex],fieldLocal)

  times = numpy.append(times,localTime)
  
if(times[-1] < 1/24.):
  timeUnit = 's'
  times *= 3600.*24.
elif(times[-1] < 1.):
  timeUnit = 'hrs'
  times *= 24.
elif(times[-1] < 365):
  timeUnit = 'days'
else:
  timeUnit = 'yrs'
  times /= 365

for varIndex in range(len(varNames)):
  plt.figure(varIndex+1)
  plt.plot(times,fields[varIndex])

for varIndex in range(len(varNames)):
  plt.figure(varIndex+1)
  plt.xlabel('time (%s)'%timeUnit)
  plt.title(varNames[varIndex])
  plt.savefig('statsPlots/%s%02i.png'%(varNames[varIndex],iteration))
  plt.close()
