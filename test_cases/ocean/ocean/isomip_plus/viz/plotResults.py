#!/usr/bin/env python
import numpy
from netCDF4 import Dataset

from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors

import os
import os.path

import copy

def makeFerretColormap():
  red = numpy.array([[0,0.6],
                     [0.15,1],
                     [0.35,1],
                     [0.65,0],
                     [0.8,0],
                     [1,0.75]])

  green = numpy.array([[0,0],
                       [0.1,0],
                       [0.35,1],
                       [1,0]])


  blue = numpy.array([[0,0],
                     [0.5,0],
                     [0.9,0.9],
                     [1,0.9]])

  colorCount = 21
  ferretColorList = numpy.ones((colorCount,4),float)
  ferretColorList[:,0] = numpy.interp(numpy.linspace(0,1,colorCount),red[:,0],red[:,1])
  ferretColorList[:,1] = numpy.interp(numpy.linspace(0,1,colorCount),green[:,0],green[:,1])
  ferretColorList[:,2] = numpy.interp(numpy.linspace(0,1,colorCount),blue[:,0],blue[:,1])
  ferretColorList = ferretColorList[::-1,:]

  cmap = colors.LinearSegmentedColormap.from_list('ferret',ferretColorList,N=255)
  return cmap


def computeCellPatches(mask):
  patches = []
  for iCell in range(nCells):
    if(not mask[iCell]):
      continue
    nVert = nVerticesOnCell[iCell]
    vertexIndices = verticesOnCell[iCell,:nVert]
    vertices = numpy.zeros((nVert,2))
    vertices[:,0] = 1e-3*xVertex[vertexIndices]
    vertices[:,1] = 1e-3*yVertex[vertexIndices]
    #middle = numpy.mean(vertices,axis=0)

    #expansion = 1.05
    # expand by 1% to avoid annoying gaps
    #for iVert in range(nVert):
      #vertices[iVert,:] = expansion*(vertices[iVert,:]-middle) + middle
    polygon = Polygon(vertices, True)
    patches.append(polygon)
    
  p = PatchCollection(patches, cmap=ferretMap, alpha=1.)
    
  return p
  
def plotHorizField(field, title, prefix, oceanDomain=True, vmin=None, vmax=None, figsize=[9,6]):
  outFileName = '%s/%s_%04i.png'%(options.outImageFolder,prefix,tIndex+1)
  if(os.path.exists(outFileName)):
    return
  if(oceanDomain):
    localPatches = copy.copy(oceanPatches)
    localPatches.set_array(field[oceanMask])
  else:
    localPatches = copy.copy(cavityPatches)
    localPatches.set_array(field[cavityMask])
  localPatches.set_edgecolor('face')
  localPatches.set_clim(vmin=vmin, vmax=vmax)
  
  plt.figure(figsize=figsize)
  ax = plt.subplot('111')
  ax.add_collection(localPatches)
  plt.colorbar(localPatches)
  plt.axis([0,500,0,1000])
  ax.set_aspect('equal')
  ax.autoscale(tight=True)
  plt.title(title)
  plt.savefig(outFileName)
  plt.close()  

def plotVertField(field, title, prefix, vmin=None, vmax=None, figsize=[9,6], inX=None, inZ=None):
  outFileName = '%s/%s_%04i.png'%(options.outImageFolder,prefix,tIndex+1)
  if(os.path.exists(outFileName)):
    return
  if(inX is None):
    inX = X
  if(inZ is None):
    inZ = Z
  plt.figure(figsize=figsize)
  ax = plt.subplot('111')
  plt.pcolor(1e-3*inX,inZ,field,vmin=vmin,vmax=vmax,cmap=ferretMap)
  plt.colorbar()
  ax.autoscale(tight=True)
  plt.ylim([numpy.amin(inZ),20])
  plt.title(title)
  plt.savefig(outFileName)
  plt.close()  

def plotHorizVertField(field, name, units, prefix, oceanDomain=True, vmin=None, vmax=None):
  if(vmin is None):
    vmin = numpy.amin(field)
  if(vmax is None):
    vmax = numpy.amax(field)

  print name, numpy.amin(field[cellMask == 1.0]), numpy.amax(field[cellMask == 1.0])
  plotHorizField(field[:,0], 'top %s (%s)'%(name,units), 'top%s'%prefix, oceanDomain, vmin=vmin, vmax=vmax)
  botField = numpy.zeros(nCells)
  for iCell in range(nCells):
    k = maxLevelCell[iCell]-1
    if (k < 0):
      continue
    botField[iCell] = field[iCell,k]
  mask = maxLevelCell > 0
  print 'bot.', name, numpy.amin(botField[mask]), numpy.amax(botField[mask])
  plotHorizField(botField, 'bot. %s (%s)'%(name,units), 'bot%s'%prefix, oceanDomain, vmin=vmin, vmax=vmax)
  field = field[sectionCellIndices,:].T
  plotVertField(field, '%s along center line (%s)'%(name,units), 'center%s'%prefix, vmin=vmin, vmax=vmax)

def computeSectionCellIndices():
  y = options.sectionY
  xMin = numpy.amin(xCell)
  xMax = numpy.amax(xCell)
  xs = numpy.linspace(xMin,xMax,10000)
  cellIndices = []
  for x in xs:
    distanceSquared = (x - xCell)**2 + (y-yCell)**2
    index = numpy.argmin(distanceSquared)
    if(len(cellIndices) == 0 or cellIndices[-1] != index):
      cellIndices.append(index)

  return numpy.array(cellIndices)

def cellToSectionEdges(field):
  nx = len(sectionCellIndices)
  fieldMid = field[sectionCellIndices]
  fieldEdge = numpy.ma.masked_all(nx+1)
  fieldEdge[1:-1] = 0.5*(fieldMid[0:-1]+fieldMid[1:])
  # extrapolate ends
  fieldEdge[0] = 2*fieldMid[0]-fieldEdge[1]
  fieldEdge[-1] = 2*fieldMid[-1]-fieldEdge[-2]
  return fieldEdge

parser = OptionParser()
           
parser.add_option("--outImageFolder", type="string", default="plots", dest="outImageFolder")
parser.add_option("--inFolder", type="string", default=".", dest="inFolder")
parser.add_option("--expt", type="int", default="1", dest="expt")
parser.add_option("--sectionY", type="float", default=40e3, dest="sectionY")

options, args = parser.parse_args()

rho_sw = 1026.
rho_fw = 1000.

ferretMap = makeFerretColormap()

try:
  os.makedirs(options.outImageFolder)
except OSError as e:
  pass

inFileName = '%s/init.nc'%(options.inFolder)
print inFileName
inFile = Dataset(inFileName,'r')
#sspRef = numpy.array(inFile.variables['seaSurfacePressure'])
oceanThickness = numpy.sum(inFile.variables['layerThickness'][0,:,:],axis=1)
bottomDepth = numpy.array(inFile.variables['bottomDepth'])
landIceFraction = inFile.variables['landIceFraction'][0,:]
sshRef = oceanThickness-bottomDepth
inFile.close()

inFileName = '%s/output.nc'%(options.inFolder)
outputFile = Dataset(inFileName,'r')
inFileName = '%s/barotropicStreamfunction.nc'%(options.inFolder)
bsfFile = Dataset(inFileName,'r')
inFileName = '%s/overturningStreamfunction.nc'%(options.inFolder)
useOSF = os.path.exists(inFileName)
if(useOSF):
  osfFile = Dataset(inFileName,'r')

inFileName = '%s/land_ice_fluxes.nc'%(options.inFolder)
landIceFluxesFile = Dataset(inFileName,'r')


xCell = numpy.array(outputFile.variables['xCell'])
yCell = numpy.array(outputFile.variables['yCell'])
bottomDepth = numpy.array(outputFile.variables['bottomDepth'])

nVertices = len(outputFile.dimensions['nVertices'])
nCells = len(outputFile.dimensions['nCells'])
nEdges = len(outputFile.dimensions['nEdges'])
nVertLevels = len(outputFile.dimensions['nVertLevels'])
nTime = len(outputFile.dimensions['Time'])
nTime = min(nTime,len(bsfFile.dimensions['Time']))
if(useOSF):
  nTime = min(nTime,len(osfFile.dimensions['Time']))

nVerticesOnCell = numpy.array(outputFile.variables['nEdgesOnCell'])
verticesOnCell = numpy.array(outputFile.variables['verticesOnCell'])-1
xVertex = numpy.array(outputFile.variables['xVertex'])
yVertex = numpy.array(outputFile.variables['yVertex'])

maxLevelCell = numpy.array(outputFile.variables['maxLevelCell'])-1
oceanMask = maxLevelCell >= 0
cavityMask = numpy.logical_and(oceanMask,landIceFraction > 0.01)
cellMask = numpy.zeros((nCells, nVertLevels))
for iCell in range(nCells):
  k = maxLevelCell[iCell]
  if (k < 0):
    continue
  cellMask[iCell,0:k+1] = 1.0

oceanPatches = computeCellPatches(oceanMask)
cavityPatches = computeCellPatches(cavityMask)

for tIndex in range(nTime):
  bsf = numpy.array(bsfFile.variables['barotropicStreamfunctionCell'])[tIndex,:]
  if(options.expt == 1):
    vmin=-1
    vmax=1
  else:
    vmin=-0.5
    vmax=0.5
  plotHorizField(bsf, 'barotropic streamfunction (Sv)', 'bsf', oceanDomain=True, vmin=vmin, vmax=vmax)
  
bsfFile.close()

if(useOSF):
  osfX = numpy.array(osfFile.variables['x'])
  osfZ = numpy.array(osfFile.variables['z'])

sectionCellIndices = computeSectionCellIndices()

x = cellToSectionEdges(xCell)
nx = len(x)
X = numpy.zeros((nVertLevels+1,nx))
for zIndex in range(nVertLevels+1):
  X[zIndex,:] = x

for tIndex in range(nTime):
  print tIndex+1, '/', nTime
  layerThickness = numpy.ma.masked_array(outputFile.variables['layerThickness'][tIndex,:,:],cellMask == 0.0)
  Z = numpy.ma.masked_all((nVertLevels+1,nx))
  Z[0,:] = 0.0
  for zIndex in range(0,nVertLevels):
    layerThicknessSection = cellToSectionEdges(layerThickness[:,zIndex])
    Z[zIndex+1,:] = Z[zIndex,:] - layerThicknessSection
  ZMin = numpy.ma.amin(Z,axis=0)
  offset = cellToSectionEdges(-bottomDepth) - ZMin
  for zIndex in range(nVertLevels+1):
    Z[zIndex,:] += offset

  try:
    secPerYear = 365*24*60*60
    freshwaterFlux = landIceFluxesFile.variables['landIceFreshwaterFlux'][tIndex,:]
    meltRate = freshwaterFlux/rho_fw*secPerYear
    plotHorizField(meltRate, 'melt rate (m/yr)', 'meltRate', oceanDomain=False, vmin=-100., vmax=100.)
  except KeyError:
    print "Key landIceFreshwaterFlux not found."
    pass

  try:
    flux = landIceFluxesFile.variables['landIceHeatFlux'][tIndex,:]
    plotHorizField(flux, 'ocean heat flux (W/s)', 'oceanHeatFlux', oceanDomain=False, vmin=-1e3, vmax=1e3)
  except KeyError:
    print "Key landIceHeatFlux not found."
    pass

  try:
    flux = landIceFluxesFile.variables['heatFluxToLandIce'][tIndex,:]
    plotHorizField(flux, 'ice heat flux (W/s)', 'iceHeatFlux', oceanDomain=False, vmin=-1e1, vmax=1e1)
  except KeyError:
    print "Key heatFluxToLandIce not found."
    pass

  try:
    Ti = landIceFluxesFile.variables['landIceInterfaceTemperature'][tIndex,:]
    To = landIceFluxesFile.variables['landIceBoundaryLayerTemperature'][tIndex,:]
    thermalDriving = To-Ti
    plotHorizField(thermalDriving, 'thermal driving (deg C)', 'thermalDriving', oceanDomain=False, vmin=-2, vmax=2)
  except KeyError:
    print "Key landIceInterfaceTemperature or landIceBoundaryLayerTemperature not found."
    pass
  try:
    Si = landIceFluxesFile.variables['landIceInterfaceSalinity'][tIndex,:]
    So = landIceFluxesFile.variables['landIceBoundaryLayerSalinity'][tIndex,:]
    halineDriving = So-Si
    plotHorizField(halineDriving, 'haline driving (PSU)', 'halineDriving', oceanDomain=False, vmin=-10, vmax=10)
  except KeyError:
    print "Key landIceInterfaceSalinity or landIceBoundaryLayerSalinity not found."
    pass
  try:
    uStar = landIceFluxesFile.variables['landIceFrictionVelocity'][tIndex,:]
    plotHorizField(uStar, 'friction velocity (m/s)', 'fricVel', oceanDomain=True, vmin=0, vmax=0.05)
  except KeyError:
    print "Key landIceFrictionVelocity not found."
    pass
  temp = outputFile.variables['temperature'][tIndex,:,:]
  plotHorizVertField(temp, 'temperature', 'deg C', 'Temp', oceanDomain=True, vmin=-2.5, vmax=1.0)

  salt = outputFile.variables['salinity'][tIndex,:,:]
  plotHorizVertField(salt, 'salinity', 'PSU', 'Salinity', oceanDomain=True, vmin=33.8, vmax=34.7)

  plotHorizVertField(layerThickness, 'layer thickness', 'm', 'LayerThickness', oceanDomain=True, vmin=0.0, vmax=25.0)

  ssh = outputFile.variables['ssh'][tIndex,:]
  delta_ssh = ssh-sshRef
  print 'delta_ssh', numpy.amin(delta_ssh), numpy.amax(delta_ssh)
  #oceanThickness = numpy.sum(layerThickness,axis=1)
  #sshDiff = ssh - (oceanThickness-bottomDepth)
  #print 'sshDiff', numpy.amin(sshDiff), numpy.amax(sshDiff)
  plotHorizField(delta_ssh, 'change in ssh (m)', 'deltaSSH', oceanDomain=True, vmin=-2, vmax=10)

  #rho = numpy.array(outputFile.variables['potentialDensity'])[tIndex,:,:]-1000.
  #plotHorizVertField(rho, 'potential density', 'kg/m^3 - 1000.', 'PotRho', vmin=27., vmax=28.)

#  N = numpy.array(outputFile.variables['BruntVaisalaFreqTop'])[tIndex,:,:]
#  plotHorizVertField(N, 'Brunt Vaisala freq.', '1/s', 'BruntVaisala')
  vx = numpy.array(outputFile.variables['velocityX'])[tIndex,:,:]
  plotHorizVertField(vx, 'X velocity', 'm/s', 'Vx', oceanDomain=True)
  vy = numpy.array(outputFile.variables['velocityY'])[tIndex,:,:]
  plotHorizVertField(vy, 'Y velocity', 'm/s', 'Vy', oceanDomain=True)
#  Ri = numpy.array(outputFile.variables['RiTopOfCell'])[tIndex,:,:-1]
#  plotHorizVertField(Ri, 'Richardson number', 'nondim.', 'Ri',vmin=-1., vmax=1.)
#  try:
#    vertVisc = numpy.array(outputFile.variables['vertViscTopOfCell'])[tIndex,:,:-1]
#    plotHorizVertField(vertVisc, 'Vertical viscosity', 'm^2/s', 'VertVisc')
#  except KeyError:
#    print "Key vertViscTopOfCell not found."
#    pass

#  try:
#    vertDiff = numpy.array(outputFile.variables['vertDiffTopOfCell'])[tIndex,:,:-1]
#    plotHorizVertField(vertDiff, 'Vertical difusivity', 'm^2/s', 'VertDiff')
#  except KeyError:
#    print "Key vertDiffTopOfCell not found."
#    pass

  if(useOSF):
    osf = osfFile.variables['overturningStreamfunction'][tIndex,:,:]
    if(options.expt == 1):
      vmin=-0.3
      vmax=0.3
    else:
      vmin=-0.3
      vmax=0.3
    plotVertField(osf, 'overturning streamfunction (Sv)', 'osf', vmin=vmin, vmax=vmax, inX=osfX, inZ=osfZ)
if(useOSF):
  osfFile.close()
outputFile.close()
landIceFluxesFile.close()
