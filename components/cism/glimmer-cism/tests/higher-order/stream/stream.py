#!/usr/bin/env python

# File to generate an input netcdf file readable by Glimmer-CISM,
# using a .mat file generated using an appropriate matlab script 
#
# Sytnax: 
# ./brpoc.py bproc.config
#
# Script will parse the grid section of the config file to produce proper output.

print " "
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print "!!! Warning: This test case is still under development !!!"
print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
print " "

from numpy import array,zeros,size
import sys
from scipy.io import loadmat
try:
  from netCDF4 import Dataset
except ImportError:
  from Scientific.IO.NetCDF import NetCDFFile as Dataset

def nc_from_config(configFilename):
  from ConfigParser import ConfigParser
  parser = ConfigParser()
  parser.read(configFilename)
  nx = int(parser.get("grid","ewn"))
  ny = int(parser.get("grid","nsn"))
  nz = int(parser.get("grid","upn"))
  dx = float(parser.get("grid","dew"))
  dy = float(parser.get("grid","dns"))
  filename = parser.get("CF input", "name")

  print "Writing to", filename
  try:
    netCDFfile = Dataset(filename,'w',format='NETCDF3_CLASSIC')
  except TypeError:
    netCDFfile = Dataset(filename,'w')
  netCDFfile.createDimension('time',1)
  netCDFfile.createDimension('x1',nx)
  netCDFfile.createDimension('y1',ny)
  netCDFfile.createDimension('level',nz)
  netCDFfile.createDimension('x0',nx-1) # staggered grid 
  netCDFfile.createDimension('y0',ny-1)
  from numpy import arange
  x = dx*arange(nx,dtype='float32')
  y = dx*arange(ny,dtype='float32')
  netCDFfile.createVariable('time','f',('time',))[:] = [0]
  netCDFfile.createVariable('x1','f',('x1',))[:] = x
  netCDFfile.createVariable('y1','f',('y1',))[:] = y
  netCDFfile.createVariable('x0','f',('x0',))[:] = dx/2 + x[:-1] # staggered grid
  netCDFfile.createVariable('y0','f',('y0',))[:] = dy/2 + y[:-1]
  return netCDFfile, (nx,ny,nz,dx,dy)

#Parse the config file to determine how to set up the netcdf file
nc, shape = nc_from_config(sys.argv[1])

topg     = nc.createVariable('topg',    'f',('time','y1','x1'))
thk      = nc.createVariable('thk',     'f',('time','y1','x1'))
usurf    = nc.createVariable('usurf',   'f',('time','y1','x1'))
#artm     = nc.createVariable('artm',    'f',('time','y1','x1'))
#acab     = nc.createVariable('acab',    'f',('time','y1','x1'))
#bheatflx = nc.createVariable('bheatflx','f',('time','y1','x1'))
beta     = nc.createVariable('beta',    'f',('time','y0','x0'))
#minTauf     = nc.createVariable('minTauf',    'f',('time','y0','x0'))

kinbcmask = nc.createVariable('kinbcmask','f',('time','y0','x0'))
uvel   = nc.createVariable('uvel',  'f',('time','level','y0','x0'))
vvel   = nc.createVariable('vvel',  'f',('time','level','y0','x0'))
#temp      = nc.createVariable('temp',  'f',('time','level','y1','x1'))
#bwat     = nc.createVariable('bwat',   'f',('time','y1','x1'))

#Read the data:
d = loadmat('stream.mat')

#Set the fields
topg[0,:,:] = array(d['topg'],dtype='float32')
thk[0,:,:] = array(d['thck'],dtype='float32')
usurf[0,:,:] = array(d['usrf'],dtype='float32')
#artm[0,:,:] = array(d['artm'],dtype='float32')
#acab[0,:,:] = array(d['acab'],dtype='float32')
#bheatflx[0,:,:] = array(d['bheatflx'],dtype='float32')
beta[0,:,:] = array(d['beta'],dtype='float32')

#minTauf[0,:,:]   = array(d['minTauf'],dtype='float32')
kinbcmask[0,:,:] = array(d['kinbcmask'],dtype='int8')
uvel[0,:,:,:] = array(d['uvel'],dtype='float32')
vvel[0,:,:,:] = array(d['vvel'],dtype='float32')
#temp[0,:,:,:] = array(d['temp'],dtype='float32')
#bwat[0,:,:] = array(d['bwat'],dtype='float32')

nc.close()
