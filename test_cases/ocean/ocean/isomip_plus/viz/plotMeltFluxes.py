#!/usr/bin/env python
import numpy
from netCDF4 import Dataset

from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

import glob

import os

parser = OptionParser()

parser.add_option("--in_file", type="string", default='land_ice_fluxes.nc', dest="in_file")
parser.add_option("--out_dir", type="string", default='meltPlots', dest="out_dir")
parser.add_option("--ssh_max", type="float", default=0., dest='ssh_max')
parser.add_option("--tavg_start", type="float", default=0.5, dest='tavg_start')
parser.add_option("--tavg_end", type="float", default=1.0, dest='tavg_end')
           
options, args = parser.parse_args()

paramPaths = args
paramName = None
outPath = options.out_dir
paramNumericalValues = []
if(len(paramPaths) == 0):
  paramPaths = ['.']
elif(len(paramPaths) > 1):
  paramName = paramPaths[0].split('_')[-2]
  paramValues = []
  for path in paramPaths:
    paramValue = path.split('_')[-1]
    paramValues.append(paramValue)
    try:
      paramNumericalValues.append(float(paramValue))
    except ValueError:
      pass

if len(paramNumericalValues) > 0:
  # make sure the parameters are sorted by parameter value
  indices = numpy.argsort(paramNumericalValues)
  paramNumericalValues = numpy.array([paramNumericalValues[i] for i in indices])
  paramValues = [paramValues[i] for i in indices]
  paramPaths = [paramPaths[i] for i in indices]

try:
  os.makedirs(outPath)
except OSError:
  pass

rho_fw = 1000.
secPerYear = 365*24*60*60

varNames = ['meanMeltRate', 'totalMeltFlux', 'thermalDriving', 'frictionVelocity']
unit = ['m/yr', 'GT/yr', 'deg C', 'cm/s']

times = []
fields = []
legends = []
maxTime = 0.

if paramName is not None:
  validParamIndices = numpy.zeros(len(paramValues),bool)

for paramIndex in range(len(paramPaths)):
  paramPath = paramPaths[paramIndex]
  inFileName = '%s/%s'%(paramPath,options.in_file)
  if(not os.path.exists(inFileName)):
    continue

  if(paramName is not None):
    legends.append('%s = %s'%(paramName, paramValues[paramIndex]))
    print legends[len(legends)-1]
    validParamIndices[paramIndex] = True

  inFile = Dataset(inFileName,'r')

  areaCell = inFile.variables['areaCell'][:]

  timesLocal = inFile.variables['daysSinceStartOfSim'][:]

  maxTime = max(maxTime,timesLocal[-1])

  nTime = len(inFile.dimensions['Time'])
  totalMeltFlux = numpy.ma.masked_all(nTime)
  meltArea = numpy.ma.masked_all(nTime)
  thermalDrivingSum = numpy.ma.masked_all(nTime)
  frictionVelocitySum = numpy.ma.masked_all(nTime)
  timeMask = numpy.zeros(nTime,bool)
  for tIndex in range(nTime):
    freshwaterFlux = inFile.variables['landIceFreshwaterFlux'][tIndex,:]
    if(numpy.all(freshwaterFlux == 0.)):
      continue
    fraction = inFile.variables['landIceFraction'][tIndex,:]
    ssh = inFile.variables['ssh'][tIndex,:]
    thermalDriving = (inFile.variables['landIceBoundaryLayerTemperature'][tIndex,:]
                   - inFile.variables['landIceInterfaceTemperature'][tIndex,:])
    frictionVelocity = inFile.variables['landIceFrictionVelocity'][tIndex,:]

    meltRate = freshwaterFlux/rho_fw*secPerYear
    sshMask = ssh < options.ssh_max
    totalMeltFlux[tIndex] = numpy.sum(sshMask*meltRate*areaCell)
    meltArea[tIndex] = numpy.sum(sshMask*fraction*areaCell)
    thermalDrivingSum[tIndex] = numpy.sum(sshMask*thermalDriving*areaCell)
    frictionVelocitySum[tIndex] = numpy.sum(sshMask*frictionVelocity*areaCell)
    timeMask[tIndex] = True

  times.append(timesLocal[timeMask])
  fields.append([totalMeltFlux[timeMask]/meltArea[timeMask], 1e-12*totalMeltFlux[timeMask],
                 thermalDrivingSum[timeMask]/meltArea[timeMask],
                 1e2*frictionVelocitySum[timeMask]/meltArea[timeMask]])

if(maxTime < 1/24.):
  timeUnit = 's'
  timeScale = 3600.*24.
elif(maxTime < 1.):
  timeUnit = 'hrs'
  timeScale = 24.
elif(maxTime < 365):
  timeUnit = 'days'
  timeScale = 1.
else:
  timeUnit = 'yrs'
  timeScale = 1./365


if len(paramNumericalValues) > 0:
  colormap = plt.get_cmap('hsv') 
  cNorm  = colors.Normalize(vmin=numpy.amin(paramNumericalValues), vmax=numpy.amax(paramNumericalValues))
  scalarMap = cm.ScalarMappable(norm=cNorm, cmap=colormap)
  plotColors = [[0.8*color for color in scalarMap.to_rgba(paramValue)] for paramValue in paramNumericalValues[validParamIndices]]

for varIndex in range(len(varNames)):
  plt.figure(1, figsize=[16,9],dpi=100)
  for paramIndex in range(len(fields)):
    if paramName is None:
      label = None
      color = 'b'
    else:
      label = legends[paramIndex]
      color = plotColors[paramIndex]
    plt.plot(timeScale*times[paramIndex],fields[paramIndex][varIndex], color=color, label=label)
  plt.xlabel('time (%s)'%timeUnit)
  plt.ylabel('%s (%s)'%(varNames[varIndex], unit[varIndex]))
  fileName = '%s/%s.png'%(outPath,varNames[varIndex])
  if paramName is not None:
    plt.legend(loc='best')
  plt.savefig(fileName)
  plt.close()

  if len(paramNumericalValues) > 0:
    print paramName, varNames[varIndex]
    plt.figure(1, figsize=[16,9],dpi=100)
    tavgValues = []
    for paramIndex in range(len(fields)):
      timesLocal = times[paramIndex]/365.
      mask = numpy.logical_and(timesLocal >= options.tavg_start,
                               timesLocal <= options.tavg_end)
      fieldValues = fields[paramIndex][varIndex][mask]
      if len(fieldValues) > 0:
         tavgValue = numpy.mean(fieldValues)
      else:
         tavgValue = numpy.nan
      print paramNumericalValues[validParamIndices][paramIndex], tavgValue
      tavgValues.append(tavgValue)
    plt.plot(paramNumericalValues[validParamIndices], tavgValues, '.-k', markersize=20)
    plt.xlabel('%s'%paramName)
    plt.ylabel('%s (%s)'%(varNames[varIndex], unit[varIndex]))
    fileName = '%s/tavg_%s.png'%(outPath,varNames[varIndex])
    plt.savefig(fileName)
    plt.close()

