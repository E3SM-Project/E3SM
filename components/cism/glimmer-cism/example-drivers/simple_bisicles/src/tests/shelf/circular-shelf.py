#!/usr/bin/env python
# This script runs a "Circular Shelf Experiment".
# Files are written in the "output" subdirectory.
# The script performs the following three steps:
# 1. Create a netCDF input file for Glimmer.
# 2. Run Glimmer, creating a netCDF output file.
# 3. Move any additional files written by Glimmer to the "scratch" subdirectory.
# Written by Glen Granzow at the University of Montana on April 9, 2010

# Slight alterations by SFP on 2-3-11 for Glimmer-CISM 2.0 relase

import sys, os, glob, shutil, numpy
from netCDF import *
from math import sqrt, exp
from optparse import OptionParser
from ConfigParser import ConfigParser

# Parse the command line arguments
parser = OptionParser()
parser.add_option('-b','--smooth-beta',dest='smooth_beta',action='store_true',help='Use a Gaussian function for beta')
parser.add_option('-d','--dirichlet-center',dest='dirichlet_center',action='store_true',help='Apply Dirichlet boundary condition at the center')
parser.add_option('-s','--sloped',dest='sloped',action='store_true',help='Use a conically topped ice thickness')
options, args = parser.parse_args()

# Check to see if a config file was specified on the command line.
# If not, circular-shelf.config is used.
if len(args) == 0:
  configfile = 'circular-shelf.config'
elif len(args) >= 1:
  configfile = args[0]
else:
  print '\nUsage:  python circular-shelf.py [FILE.CONFIG] [-b|--smooth-beta] [-d|--dirichlet-center] [-s|--sloped]\n'
  sys.exit(0)

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
dxs =  dx/2 + x[:-1] # staggered grid
dys =  dy/2 + y[:-1] # staggered grid
netCDFfile.createVariable('x0','f',('x0',))[:] = dxs.tolist()
netCDFfile.createVariable('y0','f',('y0',))[:] = dys.tolist()

# Calculate values for the required variables.
thk  = numpy.zeros([1,ny,nx],dtype='float32')
topg = numpy.empty([1,ny,nx],dtype='float32')
beta = numpy.empty([1,ny-1,nx-1],dtype='float32')
topg[0,:,:] = -2000
# Domain size
Lx = nx*dx
Ly = ny*dy

for i in range(nx):
  x = float(i)/(nx-1) - 0.5   # -1/2 < x < 1/2 
  xx = x*Lx                   # -L/2 < xx < L/2
  for j in range(ny):
    y = float(j)/(ny-1) - 0.5 # -1/2 < y < 1/2
    r = sqrt(x*x+y*y)     # radial distance from the center
    if r < 0.44:          # Inside a circle we have
      thk[0,j,i] = 1000     # constant ice thickness unless
      if options.sloped:  # command line option specifies
        thk[0,j,i] *= (1-r) # conical top
    if options.smooth_beta:  # command line option
      if i < nx-1 and j < ny-1: # beta is on the staggered grid
        yy = y*Ly               # -L/2 < yy < L/2
        beta[0,j,i] = 1.0 + 1.0e10*exp(-(xx*xx+yy*yy)/5.0e5) # Gaussian

if not options.smooth_beta:  # command line option is NOT present
  beta[0,:,:] =  0                             # beta is 0 almost everywhere
  #beta[0,:,:] =  1                             # beta is 1 almost everywhere
  beta[0,ny/2-1:ny/2+1,nx/2-1:nx/2+1] = 1.0e10 # but large in the center

# Add a single bedrock spike in the domain center, to "ground" shelf for 
# bisicles dycore
topg[0,(ny-1)/2-1:(ny-1)/2+2,(nx-1)/2-1:(nx-1)/2+2] = -800. 

# Create the required variables in the netCDF file.
netCDFfile.createVariable('thk', 'f',('time','y1','x1'))[:] = thk.tolist()
netCDFfile.createVariable('topg','f',('time','y1','x1'))[:] = topg.tolist()
netCDFfile.createVariable('beta','f',('time','y0','x0'))[:] = beta.tolist()

if options.dirichlet_center:  # command line option
  uvelbc = netCDFfile.createVariable('uvelbc','f',('time','level','y0','x0'))
  vvelbc = netCDFfile.createVariable('vvelbc','f',('time','level','y0','x0'))
  bc = numpy.empty(1,[ny-1,nx-1],dtype='float32')
# boundary condition is NaN almost everywhere
  bc[0,:,:] = float('NaN')
# boundary condition is 0 in the center
  bc[0,ny/2-1:ny/2+2,nx/2-1:nx/2+2] = 0
  for k in range(nz): # loop over levels
    uvelbc[0,k,:,:] = bc.tolist()
    vvelbc[0,k,:,:] = bc.tolist()

netCDFfile.close()

# Run Glimmer
print 'Running Glimmer/CISM'
if len(sys.argv) > 2:
   os.system('mpirun -np '+nprocs+' ./simple_glide '+configfile+'')  # support for MPI runs is here (other)
   #os.system('aprun -n'+nprocs+' ./simple_glide '+configfile+'')  # support for MPI runs is here (Jaguar)
else:
   os.system('echo '+configfile+' | ./simple_glide')

# Clean up by moving extra files written by Glimmer to the "scratch" subdirectory
# Look for files with extension "txt", "log", or "nc"
for files in glob.glob('*.txt')+glob.glob('*.log'):
# Delete any files already in scratch with these filenames 
  if files in os.listdir('scratch'):
    os.remove(os.path.join('scratch',files))
# Move the new files to scratch
  shutil.move(files,'scratch')
