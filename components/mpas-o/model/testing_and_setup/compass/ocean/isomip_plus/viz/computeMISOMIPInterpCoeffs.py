#!/usr/bin/env python
import numpy
from netCDF4 import Dataset
from shapely.geometry import Polygon, LineString

from optparse import OptionParser

from progressbar import ProgressBar, Percentage, Bar, ETA

import os.path

def getTransectWeights(outFileName, axis):

  pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=nCells).start()

  if axis == 'x':
    slicePos = xTransect
    outNOther = outNy
    outOtherAxis = y
  else:
    slicePos = yTransect
    outNOther = outNx
    outOtherAxis = x

  sliceCount = numpy.zeros(outNOther,int)
  cellIndices = []
  otherIndices = []
  weights = []
  sliceIndices = []

  for iCell in range(nCells):
    verts = verticesOnCell[iCell,0:nEdgesOnCell[iCell]]
    verts = numpy.append(verts,verts[0])
    xVert = xVertex[verts]
    yVert = yVertex[verts]

    if axis == 'x':
      sliceAxisVerts = xVert
      otherAxisVerts = yVert
    else:
      sliceAxisVerts = yVert
      otherAxisVerts = xVert



    if(numpy.amax(sliceAxisVerts) < slicePos) or (numpy.amin(sliceAxisVerts) > slicePos):
      # this polygon doesn't intersect the slice
      continue

    mpasPolygon = Polygon(zip(xVert,yVert))

    indices = numpy.nonzero(outOtherAxis < numpy.amin(otherAxisVerts))[0]
    if len(indices) == 0:
      lower = 0
    else:
      lower = indices[-1]

    indices = numpy.nonzero(outOtherAxis > numpy.amax(otherAxisVerts))[0]
    if len(indices) == 0:
      upper = outNOther
    else:
      upper = indices[0]

    for otherIndex in range(lower, upper):
      if axis == 'x':
        vertices = ((slicePos, outOtherAxis[otherIndex]),
                    (slicePos, outOtherAxis[otherIndex+1]))
      else:
        vertices = ((outOtherAxis[otherIndex], slicePos),
                    (outOtherAxis[otherIndex+1], slicePos))

      line = LineString(vertices)
      if not mpasPolygon.intersects(line):
        continue
      length = mpasPolygon.intersection(line).length
      if length == 0.0:
        continue

      cellIndices.append(iCell)
      otherIndices.append(otherIndex)
  
      weights.append(length/outDx)
  
      sliceIndex = sliceCount[otherIndex]
      sliceCount[otherIndex] += 1
      sliceIndices.append(sliceIndex)
    pbar.update(iCell+1)
  
  pbar.finish()

  cellIndices = numpy.array(cellIndices)
  otherIndices = numpy.array(otherIndices)
  sliceIndices = numpy.array(sliceIndices)
  weights = numpy.array(weights)

  # sort the intersections first by otherIndex, then by sliceIndex
  # for efficiency
  nSlices = numpy.amax(sliceCount)
  sortedIndices = numpy.zeros(0,int)
  pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=nSlices).start()
  for sliceIndex in range(nSlices):
    intersectionsInSlice = numpy.nonzero(sliceIndices == sliceIndex)[0]
    indices = numpy.argsort(otherIndices[intersectionsInSlice])
    sortedIndices = numpy.append(sortedIndices,intersectionsInSlice[indices])
    pbar.update(sliceIndex+1)
  pbar.finish()
  
  outFile = Dataset(outFileName,'w',format='NETCDF4')
  outFile.createDimension('nIntersections', len(cellIndices))
  outFile.createVariable('cellIndices','i4',('nIntersections',))
  if axis == 'x':
    outFile.createVariable('yIndices','i4',('nIntersections',))
  else:
    outFile.createVariable('xIndices','i4',('nIntersections',))
  outFile.createVariable('sliceIndices','i4',('nIntersections',))
  outFile.createVariable('mpasToMisomipWeights','f8',('nIntersections',))
  
  outVars = outFile.variables
  outVars['cellIndices'][:] = cellIndices[sortedIndices]
  if axis == 'x':
    outVars['yIndices'][:] = otherIndices[sortedIndices]
  else:
    outVars['xIndices'][:] = otherIndices[sortedIndices]

  outVars['sliceIndices'][:] = sliceIndices[sortedIndices]
  outVars['mpasToMisomipWeights'][:] = weights[sortedIndices]
  
  outFile.close()

parser = OptionParser()
options, args = parser.parse_args()

if(len(args) == 0):
  folder = '.'
else:
  folder = args[0]
meshFileName = '%s/init.nc'%folder
interpWeightsFileName = '%s/intersections.nc'%folder
xTransectFileName = '%s/xTransectIntersections.nc'%folder
yTransectFileName = '%s/yTransectIntersections.nc'%folder

outDx = 2e3
outX0 = 320e3
outY0 = 0.
xTransect = 520.01e3 # shifted slightly to avoid exact touching
yTransect = 40.01e3
outDz = 5.

outNx = 240
outNy = 40
outNz = 144

# x, y and z of corners of grid cells
x = outX0 + outDx*(numpy.arange(outNx+1))
y = outY0 + outDx*(numpy.arange(outNy+1))
z = -outDz*(numpy.arange(outNz+1))

inFile = Dataset(meshFileName,'r')

nCells = len(inFile.dimensions['nCells'])
nVertices = len(inFile.dimensions['nVertices'])
inVars = inFile.variables
nEdgesOnCell = inVars['nEdgesOnCell'][:]
verticesOnCell = inVars['verticesOnCell'][:,:]-1
xVertex = inVars['xVertex'][:]
yVertex = inVars['yVertex'][:]

inFile.close()
if(not os.path.exists(interpWeightsFileName)):
  
  cellIndices = []
  xIndices = []
  yIndices = []
  weights = []
  sliceIndices = []
  
  sliceCount = numpy.zeros((outNy,outNx),int)
  
  pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=nCells).start()
  
  for iCell in range(nCells):
    verts = verticesOnCell[iCell,0:nEdgesOnCell[iCell]]
    verts = numpy.append(verts,verts[0])
    xVert = xVertex[verts]
    yVert = yVertex[verts]
    mpasPolygon = Polygon(zip(xVert,yVert))
    #mpasArea = mpasPolygon.area
  
    # find the out indices that bound the MPAS polygon
    indices = numpy.nonzero(x < numpy.amin(xVert))[0]
    if len(indices) == 0:
      xl = 0
    else:
      xl = indices[-1]
  
    indices = numpy.nonzero(x > numpy.amax(xVert))[0]
    if len(indices) == 0:
      xu = outNx
    else:
      xu = indices[0]
  
    indices = numpy.nonzero(y < numpy.amin(yVert))[0]
    if len(indices) == 0:
      yl = 0
    else:
      yl = indices[-1]
  
    indices = numpy.nonzero(y > numpy.amax(yVert))[0]
    if len(indices) == 0:
      yu = outNy
    else:
      yu = indices[0]
  
    for yIndex in range(yl, yu):
      for xIndex in range(xl, xu):
        vertices = ((x[xIndex], y[yIndex]),
                    (x[xIndex+1], y[yIndex]),
                    (x[xIndex+1], y[yIndex+1]),
                    (x[xIndex], y[yIndex+1]),
                    (x[xIndex], y[yIndex]))
        outPoly = Polygon(vertices)
        if not mpasPolygon.intersects(outPoly):
          continue
        intersectionArea = mpasPolygon.intersection(outPoly).area
        if intersectionArea == 0.:
          continue
  
        cellIndices.append(iCell)
        xIndices.append(xIndex)
        yIndices.append(yIndex)
  
        weights.append(intersectionArea/outDx**2)
  
        sliceIndex = sliceCount[yIndex,xIndex]
        sliceCount[yIndex,xIndex] += 1
        sliceIndices.append(sliceIndex)
    pbar.update(iCell+1)
  
  pbar.finish()
  
  cellIndices = numpy.array(cellIndices)
  xIndices = numpy.array(xIndices)
  yIndices = numpy.array(yIndices)
  sliceIndices = numpy.array(sliceIndices)
  weights = numpy.array(weights)
  
  # sort the intersections first by xIndex, then by yIndex, then by sliceIndex
  # for efficiency
  nSlices = numpy.amax(sliceCount)
  sortedIndices = numpy.zeros(0,int)
  xyIndices = xIndices + outNx*yIndices
  pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=nSlices).start()
  for sliceIndex in range(nSlices):
    intersectionsInSlice = numpy.nonzero(sliceIndices == sliceIndex)[0]
    indices = numpy.argsort(xyIndices[intersectionsInSlice])
    sortedIndices = numpy.append(sortedIndices,intersectionsInSlice[indices])
    pbar.update(sliceIndex+1)
  pbar.finish()
  
  outFile = Dataset(interpWeightsFileName,'w',format='NETCDF4')
  outFile.createDimension('nIntersections', len(cellIndices))
  outFile.createVariable('cellIndices','i4',('nIntersections',))
  outFile.createVariable('xIndices','i4',('nIntersections',))
  outFile.createVariable('yIndices','i4',('nIntersections',))
  outFile.createVariable('sliceIndices','i4',('nIntersections',))
  outFile.createVariable('mpasToMisomipWeights','f8',('nIntersections',))
  
  outVars = outFile.variables
  outVars['cellIndices'][:] = cellIndices[sortedIndices]
  outVars['xIndices'][:] = xIndices[sortedIndices]
  outVars['yIndices'][:] = yIndices[sortedIndices]
  outVars['sliceIndices'][:] = sliceIndices[sortedIndices]
  outVars['mpasToMisomipWeights'][:] = weights[sortedIndices]
  
  outFile.close()

if(not os.path.exists(xTransectFileName)):
  getTransectWeights(xTransectFileName, axis='x')

if(not os.path.exists(yTransectFileName)):
  getTransectWeights(yTransectFileName, axis='y')
  
