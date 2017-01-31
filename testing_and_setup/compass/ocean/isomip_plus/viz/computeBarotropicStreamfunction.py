#!/usr/bin/env python
import numpy
from netCDF4 import Dataset

from optparse import OptionParser
#import matplotlib.pyplot as plt

import scipy.sparse
import scipy.sparse.linalg
import os.path


def computeTransport():
  transport = numpy.zeros(nEdges)
  innerEdges = numpy.logical_and(cellsOnEdge[:,0] >= 0,cellsOnEdge[:,1] >= 0)
  innerEdges = numpy.nonzero(innerEdges)[0]
  for iEdge in innerEdges:
    cell0 = cellsOnEdge[iEdge,0]
    cell1 = cellsOnEdge[iEdge,1]
    layerThicknessEdge = 0.5*(layerThickness[cell0,:] + layerThickness[cell1,:])
    transport[iEdge] = dvEdge[iEdge]*numpy.sum(layerThicknessEdge*normalVelocity[iEdge,:])

  return (innerEdges,transport)
  
def computeBSF():
  boundaryVertices = numpy.logical_or(cellsOnVertex[:,0] ==-1,cellsOnVertex[:,1] ==-1)
  boundaryVertices = numpy.logical_or(boundaryVertices,cellsOnVertex[:,2] ==-1)
  boundaryVertices = numpy.nonzero(boundaryVertices)[0]
  nBoundaryVertices = len(boundaryVertices)
  nInnerEdges = len(innerEdges)
  
  rhs = numpy.zeros(nInnerEdges+nBoundaryVertices,dtype=float)
  indices = numpy.zeros((2,2*nInnerEdges+nBoundaryVertices),dtype=int)
  data = numpy.zeros(2*nInnerEdges+nBoundaryVertices,dtype=float)
  for index in range(nInnerEdges):
  
    iEdge = innerEdges[index]
    v0 = verticesOnEdge[iEdge,0]
    v1 = verticesOnEdge[iEdge,1]
    
    indices[0,2*index] = index
    indices[1,2*index] = v1
    data[2*index] = 1.
    indices[0,2*index+1] = index
    indices[1,2*index+1] = v0
    data[2*index+1] = -1.
    
    rhs[index] = transport[iEdge]*1e-6 #in Sv
  
  for index in range(nBoundaryVertices):
    iVertex = boundaryVertices[index]
    indices[0,2*nInnerEdges+index] = nInnerEdges+index
    indices[1,2*nInnerEdges+index] = iVertex
    data[2*nInnerEdges+index] = 1. 
    rhs[nInnerEdges+index] = 0. # bsf is zero at the boundaries
  
  M = scipy.sparse.csr_matrix((data,indices),shape=(nInnerEdges+nBoundaryVertices,nVertices))
  
  solution = scipy.sparse.linalg.lsqr(M,rhs)
  
  bsf = -solution[0]
  return bsf

 
def computeBSFCell(bsf):
  bsfCell = numpy.zeros(nCells)
  for iCell in range(nCells):
    edgeCount = nEdgesOnCell[iCell]
    edges = edgesOnCell[iCell,0:edgeCount]
    verts = verticesOnCell[iCell,0:edgeCount]
    areaEdge = dcEdge[edges]*dvEdge[edges]
    indexM1 = numpy.mod(numpy.arange(-1,edgeCount-1),edgeCount)
    areaVert = 0.5*(areaEdge+areaEdge[indexM1])
    bsfCell[iCell] = numpy.sum(areaVert*bsf[verts])/numpy.sum(areaVert)
  
  return bsfCell

parser = OptionParser()
           
options, args = parser.parse_args()

folder=args[0]

inFile = Dataset('%s/output.nc'%folder,'r')
outFileName = '%s/barotropicStreamfunction.nc'%folder
continueOutput = os.path.exists(outFileName)
if(continueOutput):
  outFile = Dataset(outFileName,'r+',format='NETCDF4')
else:
  outFile = Dataset(outFileName,'w',format='NETCDF4')

  #Copy dimensions
  for dimName, dimension in inFile.dimensions.iteritems():
    print dimName, len(dimension)
    outFile.createDimension(dimName, len(dimension) if not dimension.isunlimited() else None)
 
 
  # Copy variables
  for varName, inVar in inFile.variables.iteritems():
    if(inVar.dimensions[0] == 'Time'):
      continue
    outVar = outFile.createVariable(varName, inVar.datatype, inVar.dimensions)
    print varName, inVar.dimensions, inVar.datatype
  
    # Copy variable attributes
    outVar.setncatts({k: inVar.getncattr(k) for k in inVar.ncattrs()})
  
    outVar[:] = inVar[:]

nVertices = len(inFile.dimensions['nVertices'])
nCells = len(inFile.dimensions['nCells'])
nEdges = len(inFile.dimensions['nEdges'])
nVertLevels = len(inFile.dimensions['nVertLevels'])
nTimeIn = len(inFile.dimensions['Time'])
if(continueOutput):
  nTimeOut = len(outFile.dimensions['Time'])
else:
  nTimeOut = 0

verticesOnEdge = inFile.variables['verticesOnEdge'][:]-1 # change to zero-based indexing
cellsOnEdge = inFile.variables['cellsOnEdge'][:]-1 # change to zero-based indexing
dvEdge = inFile.variables['dvEdge'][:]
dcEdge = inFile.variables['dcEdge'][:]
cellsOnVertex = inFile.variables['cellsOnVertex'][:]-1 # change to zero-based indexing
nEdgesOnCell = inFile.variables['nEdgesOnCell'][:]
edgesOnCell = inFile.variables['edgesOnCell'][:]-1
verticesOnCell = inFile.variables['verticesOnCell'][:]-1
xEdge = inFile.variables['xEdge'][:]
yEdge = inFile.variables['yEdge'][:]
xVertex = inFile.variables['xVertex'][:]
yVertex = inFile.variables['yVertex'][:]

if(continueOutput):
  outBSF = outFile.variables['barotropicStreamfunction']
  outBSFCell = outFile.variables['barotropicStreamfunctionCell']
else:
  outBSF = outFile.createVariable('barotropicStreamfunction',float,['Time','nVertices'])
  outBSFCell = outFile.createVariable('barotropicStreamfunctionCell',float,['Time','nCells'])

print nTimeOut, nTimeIn
for tIndex in range(nTimeOut,nTimeIn):
  print tIndex, nTimeIn
  normalVelocity = inFile.variables['normalVelocity'][tIndex,:,:]
  layerThickness = inFile.variables['layerThickness'][tIndex,:,:]
  

  (innerEdges,transport) = computeTransport()
  
  bsf = computeBSF()
  
  bsfCell = computeBSFCell(bsf)
  
  
  outBSF[tIndex,:] = bsf
  outBSFCell[tIndex,:] = bsfCell

outFile.close()
inFile.close()

