#!/usr/bin/env python
import numpy
from netCDF4 import Dataset

from optparse import OptionParser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

import os
import os.path

import copy

def computeCellPatches():
  patches = []
  for iCell in range(nCells):
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

  p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=1.)

  return p

def plotHorizField(field, title, prefix, vmin=None, vmax=None, figsize=[6,9]):
  outFileName = '%s/%s_%04i.png'%(options.outImageFolder,prefix,tIndex+1)
  #if(os.path.exists(outFileName)):
  #  return

  if(vmin is None):
    vmin = numpy.amin(field)
  if(vmax is None):
    vmax = numpy.amax(field)
  localPatches = copy.copy(cellPatches)
  localPatches.set_array(field)
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

def plotVertField(field, title, prefix, vmin=None, vmax=None, figsize=[9,6], inY=None, inZ=None):
  outFileName = '%s/%s_%04i.png'%(options.outImageFolder,prefix,tIndex+1)
  if(os.path.exists(outFileName)):
    return
  if(inY is None):
    inY = Y
  if(inZ is None):
    inZ = Z
  plt.figure(figsize=figsize)
  ax = plt.subplot('111')
  plt.pcolor(1e-3*inY,inZ,field,vmin=vmin,vmax=vmax)
  plt.colorbar()
  ax.autoscale(tight=True)
  plt.ylim([numpy.amin(inZ),0])
  plt.title(title)
  plt.savefig(outFileName)
  plt.close()

def plotHorizVertField(field, name, units, prefix, vmin=None, vmax=None):
  if(vmin is None):
    vmin = numpy.amin(field)
  if(vmax is None):
    vmax = numpy.amax(field)

  print name, numpy.amin(field), numpy.amax(field)
  plotHorizField(field[:,0], 'top %s (%s)'%(name,units), 'top%s'%prefix, vmin=vmin, vmax=vmax)
  field = field[sectionCellIndices,:].T
  plotVertField(field, '%s along center line (%s)'%(name,units), 'center%s'%prefix, vmin=vmin, vmax=vmax)

def computeSectionCellIndices():
  x = options.sectionX
  yMin = numpy.amin(yCell)
  yMax = numpy.amax(yCell)
  ys = numpy.linspace(yMin,yMax,10000)
  cellIndices = []
  for y in ys:
    distanceSquared = (x - xCell)**2 + (y-yCell)**2
    index = numpy.argmin(distanceSquared)
    if(len(cellIndices) == 0 or cellIndices[-1] != index):
      cellIndices.append(index)

  return numpy.array(cellIndices)

def cellToSectionEdges(field):
  ny = len(sectionCellIndices)
  fieldMid = field[sectionCellIndices]
  fieldEdge = numpy.zeros(ny+1)
  fieldEdge[1:-1] = 0.5*(fieldMid[0:-1]+fieldMid[1:])
  # extrapolate ends
  fieldEdge[0] = 2*fieldMid[0]-fieldEdge[1]
  fieldEdge[-1] = 2*fieldMid[-1]-fieldEdge[-2]
  return fieldEdge

parser = OptionParser()

parser.add_option("--outImageFolder", type="string", default="plots", dest="outImageFolder")
parser.add_option("--inFolder", type="string", default=".", dest="inFolder")
parser.add_option("--initFile", type="string", default='init.nc', dest="initFile")
parser.add_option("--sshFile", type="string", default='output_ssh.nc', dest="sshFile")
parser.add_option("--iterIndex", type="int", default=0, dest="iterIndex")

options, args = parser.parse_args()

try:
  os.makedirs(options.outImageFolder)
except OSError as e:
  pass

inFileName = '%s/%s'%(options.inFolder,options.initFile)
inFile = Dataset(inFileName,'r')
inFileName = '%s/%s'%(options.inFolder,options.sshFile)
sshFile = Dataset(inFileName,'r')

nVertices = len(inFile.dimensions['nVertices'])
nCells = len(inFile.dimensions['nCells'])
nEdges = len(inFile.dimensions['nEdges'])
nVertLevels = len(inFile.dimensions['nVertLevels'])

nVerticesOnCell = numpy.array(inFile.variables['nEdgesOnCell'])
verticesOnCell = numpy.array(inFile.variables['verticesOnCell'])-1
xVertex = numpy.array(inFile.variables['xVertex'])
yVertex = numpy.array(inFile.variables['yVertex'])

nTime = len(sshFile.dimensions['Time'])
landIcePressure = inFile.variables['landIcePressure'][nTime-1,:]
nTime = len(sshFile.dimensions['Time'])
ssh = sshFile.variables['ssh'][nTime-1,:]
deltaSSH = ssh - inFile.variables['ssh'][0,:]

inFile.close()

cellPatches = computeCellPatches()

tIndex = options.iterIndex

plotHorizField(landIcePressure, 'land-ice pressure (Pa)', 'landIcePressure')
plotHorizField(ssh, 'SSH (m)', 'ssh')
plotHorizField(deltaSSH, 'delta SSH (m)', 'deltaSSH')

