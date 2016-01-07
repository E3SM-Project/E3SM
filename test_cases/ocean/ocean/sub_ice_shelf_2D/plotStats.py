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

varNames = ['kineticEnergyCellMax', 'kineticEnergyCellAvg','layerThicknessMin']
timeUnit = 'hours'
sPerUnit = 60*60

field1 = None

firstYear = None
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

  nTime = len(inFile.dimensions['Time'])
  for varIndex in range(len(varNames)):
    fieldLocal = numpy.array(inFile.variables[varNames[varIndex]])
    fields[varIndex] = numpy.append(fields[varIndex],fieldLocal)
  xTimes = numpy.array(inFile.variables['xtime'])
  localTime = numpy.zeros(nTime)
  for tIndex in range(len(xTimes)):
    xTime = ''.join(xTimes[tIndex]).strip()
    year = int(xTime[0:4])
    if firstYear is None:
      firstYear = year
      offsetYear = (firstYear == 0)
      if(offsetYear):
        firstYear += 1
      firstTime = time.mktime(time.strptime('%04i'%firstYear, '%Y'))
    
    if(offsetYear):
      xTime = '%04i'%(year+1)+xTime[4:]

    timeStruct = time.strptime(xTime, '%Y-%m-%d_%H:%M:%S')
    localTime[tIndex] = (time.mktime(timeStruct)-firstTime)/sPerUnit
    #print localTime[tIndex]
  times = numpy.append(times,localTime)
  
for varIndex in range(len(varNames)):
  plt.figure(varIndex+1)
  plt.plot(times,fields[varIndex])

for varIndex in range(len(varNames)):
  plt.figure(varIndex+1)
  plt.xlabel('time (%s)'%timeUnit)
  plt.title(varNames[varIndex])
  plt.savefig('statsPlots/%s%02i.png'%(varNames[varIndex],iteration))
  plt.close()
