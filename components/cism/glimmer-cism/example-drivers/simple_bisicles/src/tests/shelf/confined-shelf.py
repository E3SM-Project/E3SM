#!/usr/bin/env python
# This script runs a "Confined Shelf Experiment".
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Written by Glen Granzow at the University of Montana on April 8, 2010

# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 relase

import sys, os, glob, shutil, numpy
from netCDF import *
from math import sqrt
from ConfigParser import ConfigParser

# Check to see if a config file was specified on the command line.
# If not, confined-shelf.config is used.
if len(sys.argv) > 1:
  if sys.argv[1][0] == '-': # The filename can't begin with a hyphen
    print '\nUsage:  python confined-shelf.py [FILE.CONFIG]\n'
    sys.exit(0)
  else:
    configfile = sys.argv[1]
else:
  configfile = 'confined-shelf.config'

# Check to see if #procs specified, relevant when running the code in parallel.
# If not, serial run (#procs==1) is performed. To run in parallel, the configure
# file must be specifed, but the nu,ber of processors does not
if len(sys.argv) > 2:
    nprocs = sys.argv[2]
else:
  nprocs = '1'

# Create a netCDF file according to the information in the config file.
parser = ConfigParser()
parser.read(configfile)
nx = int(parser.get('grid','ewn'))
ny = int(parser.get('grid','nsn'))
nz = int(parser.get('grid','upn'))
dx = float(parser.get('grid','dew'))
dy = float(parser.get('grid','dns'))
filename = parser.get('CF input', 'name')

# Domain size
Lx = nx*dx
Ly = ny*dy

print 'Writing', filename
try:
  netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
except TypeError:
  netCDFfile = NetCDFFile(filename,'w')

netCDFfile.createDimension('time',1)
netCDFfile.createDimension('x1',nx)
netCDFfile.createDimension('y1',ny)
netCDFfile.createDimension('level',nz)
netCDFfile.createDimension('x0',nx-1) # staggered grid 
netCDFfile.createDimension('y0',ny-1)

x = dx*numpy.arange(nx,dtype='float32')
y = dx*numpy.arange(ny,dtype='float32')

netCDFfile.createVariable('time','f',('time',))[:] = [0]
netCDFfile.createVariable('x1','f',('x1',))[:] = x.tolist()
netCDFfile.createVariable('y1','f',('y1',))[:] = y.tolist()

netCDFfile.createVariable('x0','f',('x0',))[:] = (dx/2 + x[:-1]).tolist()
netCDFfile.createVariable('y0','f',('y0',))[:] = (dy/2 + y[:-1]).tolist()

# *SFP* this has been changed so that the default value for 'flwa' is the same 
# as in the EISMINT-shelf test documentation, tests 3 & 4, found at:
# http://homepages.vub.ac.be/~phuybrec/eismint/iceshelf.html

## Check to make sure that the flow law parameter in the config file is correct.
#default_flwa = float(parser.get('parameters','default_flwa'))
#if default_flwa != 4.6e-18:
#  print 'WARNING: The parameter default_flwa in',configfile,'should be 4.6e-18'
#  print '         Currently it is',default_flwa

# *SFP* removed periodic option
# Determine from the config file whether periodic boundary conditions are to be
# imposed in the x direction.  
#periodic_ew = int(parser.get('options','periodic_ew'))

# Calculate values for the required variables.
thk  = numpy.zeros([1,ny,nx],dtype='float32')
topg  = numpy.zeros([1,ny,nx],dtype='float32')
beta = numpy.empty([1,ny-1,nx-1],dtype='float32')
kbc  = numpy.zeros([1,ny-1,nx-1],dtype='int')
acab = numpy.zeros([1,ny,nx],dtype='float32') # *sfp* added acab field for prog. runs 
bmlt = numpy.zeros([1,ny,nx],dtype='float32') # *dfm* added bmlt field for testing subshelf melt coupling 
temp = numpy.zeros([1,nz,ny,nx],dtype='float32') 
zero = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

uvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')
vvel = numpy.zeros([1,nz,ny-1,nx-1],dtype='float32')

# *SFP* added topg var so that areas of no slip are consistent w/ grounded ice on bedrock
topg[:] = -2000.0

# *SFP* changed to be in line w/ EISMINT-shelf tests 3&4 

# shelf bc applied at bottom (DEFAULT FOR TEST CASE - other options below for testing bcs)
thk[0,4:-2,2:-2] = 500.     
kbc[0,ny-4:,:]  = 1
kbc[0,:,:3] = 1
kbc[0,:,nx-4:] = 1
topg[0,ny-4:,:]  = -440 
topg[0,:,:4] = -440
topg[0,:,nx-4:] = -440

# shelf bc applied at top    
#thk[0,2:-4,2:-2] = 500.     
#kbc[0,:3,:]  = 1
#kbc[0,:,:3] = 1
#kbc[0,:,nx-4:] = 1
#topg[0,:4,:]  = -440
#topg[0,:,:4] = -440
#topg[0,:,nx-4:] = -440

# shelf bc applied at right     ! NOTE that shelf is wider slightly wider in ns than in ew direction  
#thk[0,2:-2,2:-4] = 500.     
#kbc[0,:,:3]  = 1
#kbc[0,:3,:] = 1
#kbc[0,ny-4:,:] = 1
#topg[0,:,:4]  = -440
#topg[0,:4,:] = -440
#topg[0,ny-4:,:] = -440

# shelf bc applied at left     ! NOTE that shelf is wider slightly wider in ns than in ew direction  
#thk[0,2:-2,4:-2] = 500.     
#kbc[0,:,nx-4:]  = 1
#kbc[0,:3,:] = 1
#kbc[0,ny-4:,:] = 1
#topg[0,:,nx-4:]  = -440
#topg[0,:4,:] = -440
#topg[0,ny-4:,:] = -440

#if not periodic_ew:    *SFP* removed periodic option

beta[0,:,:] = 0 

acab[:] = 0.25

for i in range(nx):
  x = float(i)/(nx-1) - 0.65   # -1/2 < x < 1/2 
  xx = x*Lx                   # -L/2 < xx < L/2
  for j in range(ny):
    y = float(j)/(ny-1) - 0.5 # -1/2 < y < 1/2
    r = sqrt(x*x+y*y)     # radial distance from the center
    if r < 0.15:          # Inside a circle we have
      acab[0,j,i] = -150.0     # really strong melting

acab[0,ny-3:,:]  = 0    # zero out accum at edges to avoid buildup where u=0
acab[0,:,:3] = 0
acab[0,:,nx-3:] = 0

#bmelt

bmlt[:] = 0.25

# Domain size
Lx = nx*dx
Ly = ny*dy

for i in range(nx):
  x = float(i)/(nx-1) - 0.65   # -1/2 < x < 1/2 
  xx = x*Lx                   # -L/2 < xx < L/2
  for j in range(ny):
    y = float(j)/(ny-1) - 0.5 # -1/2 < y < 1/2
    r = sqrt(x*x+y*y)     # radial distance from the center
    if r < 0.15:          # Inside a circle we have
      bmlt[0,j,i] = -25.0     # really strong melting

bmlt[0,ny-3:,:]  = 0    # zero out accum at edges to avoid buildup where u=0
bmlt[0,:,:3] = 0
bmlt[0,:,nx-3:] = 0

temp[:] = -10.0        

# *SFP* calculate stream profile for upstream end
for i in range(nx-2):
  x = float( i ) / (nx-2) - 0.5
  vvel[0,:,ny-4,i] = -1.5e3 * 1/(2*3.141592654*0.125) * numpy.exp( -x**2 / (2*0.125**2) )

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk',      'f',('time','y1','x1'))[:] = thk.tolist()
netCDFfile.createVariable('acab',     'f',('time','y1','x1'))[:] = acab.tolist()
netCDFfile.createVariable('bmlt',     'f',('time','y1','x1'))[:] = bmlt.tolist()
netCDFfile.createVariable('temp',     'f',('time','level','y1','x1'))[:] = temp.tolist()
netCDFfile.createVariable('kinbcmask','i',('time','y0','x0'))[:] = kbc.tolist()
netCDFfile.createVariable('topg',     'f',('time','y1','x1'))[:] = topg.tolist()
netCDFfile.createVariable('beta',     'f',('time','y0','x0'))[:] = beta.tolist()
netCDFfile.createVariable('uvel',  'f',('time','level','y0','x0'))[:] = zero.tolist()

# *sfp* first option below adds ice stream vel profile for kin bc at upstream end
# *sfp* ... comment out for standard test case
#netCDFfile.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = vvel.tolist()
netCDFfile.createVariable('vvel',  'f',('time','level','y0','x0'))[:] = zero.tolist()

netCDFfile.close()

# Clean up by moving extra files written by Glimmer to the "scratch" subdirectory
# Look for files with extension "txt", "log", or "nc"
for files in glob.glob('*.txt')+glob.glob('*.log'):
# Delete any files already in scratch with these filenames 
  if files in os.listdir('scratch'):
    os.remove(os.path.join('scratch',files))
# Move the new files to scratch
  shutil.move(files,'scratch')
