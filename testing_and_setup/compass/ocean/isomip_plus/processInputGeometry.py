#!/usr/bin/env python
import numpy
from netCDF4 import Dataset
import argparse
import os

def smoothGeometry(landFraction, floatingFraction, bed, draft, filterSigma):
  import scipy.ndimage.filters as filters

  # Smoothing is performed using only the topography in the portion of the grid that is ocean.
  # This prevents the kink in the ice draft across the grounding line or regions of bare bedrock
  # from influencing the smoothed topography.  (Parts of the Ross ice shelf near the Trans-Antarctic
  # Mountains are particularly troublesome if topogrpahy is smoothed across the grounding line.)

  # Unlike in POP, the calving front is smoothed as well because MPAS-O does not support a
  # sheer calving face

  threshold = 0.01 # we won't normalize bed topography or ice draft where the mask is below this threshold

  oceanFraction = 1. - landFraction
  smoothedMask = filters.gaussian_filter(oceanFraction,filterSigma,mode='constant',cval=0.)
  mask = smoothedMask > threshold

  draft = filters.gaussian_filter(draft*oceanFraction,filterSigma,mode='constant',cval=0.)
  draft[mask] /= smoothedMask[mask]
  bed = filters.gaussian_filter(bed*oceanFraction,filterSigma,mode='constant',cval=0.)
  bed[mask] /= smoothedMask[mask]

  smoothedDraftMask = filters.gaussian_filter(floatingFraction,filterSigma,mode='constant',cval=0.)
  smoothedDraftMask[mask] /= smoothedMask[mask]

  return (bed,draft,smoothedDraftMask)

def readVar(varName, defaultValue=0.0):
  field = defaultValue*numpy.ones((ny,nx),float)
  field[buffer:-buffer,buffer:-buffer] = numpy.array(inFile.variables[varName])[:,minIndex:]
  return field

def writeVar(outVarName, inVarName, field):
  outVar = outFile.createVariable(outVarName,'f8',('y','x'))
  inVar = inFile.variables[inVarName]
  outVar[:,:] = field
  outVar.setncatts({k: inVar.getncattr(k) for k in inVar.ncattrs()})


parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-i', dest='inFileName', required=True,
                    help='ISOMIP+ input geometry file')
parser.add_argument('-o', dest='outFileName', required=True,
                    help='output geometry file after smoothing, etc.')
parser.add_argument('-s', dest='filterSigma', type=float, required=True,
                    help='no. of grid cells over which to smooth the geometry')
parser.add_argument('-m', dest='minIceThickness', type=float, required=True,
                    help='minimum ice thickness, below which ice is "calved"')
parser.add_argument('--scale', dest='scale', type=float, default=1.0,
                    help='a scale factor by which to multiple the ice draft')

args = parser.parse_args()

buffer = 1

xMax = 800e3 #km
x0 = 320e3 #km

densityRatio = 918./1028.
  
inFile = Dataset(args.inFileName,'r')
x = numpy.array(inFile.variables['x'])[:]
y = numpy.array(inFile.variables['y'])[:]

deltaX = x[1]-x[0]
deltaY = y[1]-y[0]

minIndex = numpy.nonzero(x >= x0)[0][0]

nx = len(x)-minIndex+2*buffer
ny = len(y)+2*buffer

outX = numpy.array(nx)
outY = numpy.array(ny)

outX = x[minIndex] + deltaX*(-buffer + numpy.arange(nx))
outY = y[0] + deltaY*(-buffer + numpy.arange(ny))

surf = readVar('upperSurface')
draft = readVar('lowerSurface')
bed = readVar('bedrockTopography')
floatingMask = readVar('floatingMask')
groundedMask = readVar('groundedMask', defaultValue=1.0)
openOceanMask = readVar('openOceanMask')

iceThickness = surf-draft

draft *= args.scale

# take care of calving criterion
mask = numpy.logical_and(floatingMask > 0.1, 
                         iceThickness < args.minIceThickness)
surf[mask] = 0.
draft[mask] = 0.
floatingMask[mask] = 0.
openOceanMask[mask] = 1. - groundedMask[mask]

(bed,draft,smoothedDraftMask) = smoothGeometry(groundedMask, floatingMask, bed, draft, args.filterSigma)

outFile = Dataset(args.outFileName,'w', format='NETCDF4')
outFile.createDimension('x', nx)
outFile.createDimension('y', ny)

outVar = outFile.createVariable('x','f8',('x'))
inVar = inFile.variables['x']
outVar[:] = outX
outVar.setncatts({k: inVar.getncattr(k) for k in inVar.ncattrs()})
outVar = outFile.createVariable('y','f8',('y'))
inVar = inFile.variables['y']
outVar[:] = outY
outVar.setncatts({k: inVar.getncattr(k) for k in inVar.ncattrs()})
writeVar('Z_ice_surface', 'upperSurface', surf)
writeVar('Z_ice_draft', 'lowerSurface', draft)
writeVar('Z_bed', 'bedrockTopography', bed)
writeVar('floatingIceFraction', 'floatingMask', floatingMask)
writeVar('landFraction', 'groundedMask', groundedMask)
writeVar('openOceanFraction', 'openOceanMask', openOceanMask)
writeVar('smoothedDraftMask', 'openOceanMask', smoothedDraftMask)

outFile.close()
inFile.close()

