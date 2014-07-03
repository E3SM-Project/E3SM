#!/usr/bin/env python
# This script runs a Ross Ice Shelf experiment using Glimmer.
# Output files are written in the "output" subdirectory.
# Two netCDF files are created; one contains "raw" data, the other is an input file for Glimmer.
# After creating the neccessary netCDF input file, Glimmer is run.
# Glimmer writes an additional netCDF output file in the "output" subdirectory.
# When finished, any additional files written by Glimmer are moved to the "scratch" subdirectory.
# After running this script, run plotRoss.py to plot the results.
# See the accompanying README file for more information.
# Written March 18, 2010 by Glen Granzow at the University of Montana.

import os
import glob
import shutil
import numpy
import ConfigParser
from netCDF import *

create_files = True 
run_glimmer  = False
executable   = 'simple_glide'
fake_shelf   = False
verbose      = False
use_inlets   = (False, True, 'reverse')[1]
mask_llc     = False # The previous script (makerossnc.py) masked the lower left corner
offset_error = 1     # offset_error should be 1 

def addBorder(data,dtype,value=0):
  field = numpy.empty((ny,nx),dtype)
  field[:] = value
  field[2:-2,2:-2] = data
  return field

def plot(variable): # used for debugging only
  from matplotlib import pyplot
  pyplot.figure()
  pyplot.imshow(variable,origin='lower',interpolation='nearest')
  pyplot.colorbar()
  pyplot.show()

########## PART I: READ THE INPUT FILES ##########

if create_files:

# Read the main data file into a dictionary mapping names to lists (of lists)
# The dictionary keys are the headers that begin each section in the file
  filename = os.path.join('data','111by147Grid.dat')
  inputfile = open(filename)
  print '\nReading',filename
  currentKey = None
  currentList = list()
  data = dict()
  for line in inputfile:
    line = line.strip()
    if line.startswith('#'):
      if currentKey != None:
        data[currentKey] = currentList
        currentList = list()
      currentKey = line[1:].strip().lower()
    elif len(line) > 0:
      if line.find('.') > 0:
        currentList.append([float(x) for x in line.split()])
      else:
        currentList.append([int(x) for x in line.split()])
  data[currentKey] = currentList
  inputfile.close()

  print 'The',len(data.keys()),'data fields read from 111by147Grid.dat are:\n',data.keys()

# Read the kinematic boundary conditions (kbc) mask file.  
# This is a set of (i,j) coordinates specifying where the velocity read from
# the data file should be used as a boundary condition for the model.
  kbc_mask = numpy.zeros((111,147), dtype='i')
  filename = os.path.join('data','kbc.dat')
  inputfile = open(filename)
  print '\nReading',filename
  for line in inputfile:
    i,j = map(int,line.split())
    if offset_error != 0: (i,j) = (i-offset_error,j-offset_error)
    kbc_mask[i,j] = 1
  inputfile.close()
  print numpy.sum(kbc_mask),'points were read from kbc.dat'

# Read in the inlets file, which specifies additional Dirichlet boundary conditions
  filename = os.path.join('data','inlets.dat')
  inputfile = open(filename)
  print '\nReading',filename
  counter = [0,0]
  for line in inputfile:
    i, j, azimuth, magnitude = line.split()
    i,j = map(int,(i,j))
    if offset_error != 0: (i,j) = (i-offset_error,j-offset_error)
    if use_inlets == 'reverse':
      magnitude,azimuth = map(float,(azimuth,magnitude))
    else:
      azimuth,magnitude = map(float,(azimuth,magnitude))
    if use_inlets:
      if verbose:
        indices = '(%d,%d):' % (i,j)
        print 'Changing azimuth  at',indices,data['ice velocity azimuth grid'][i][j],'->',azimuth
        print 'Changing velocity at',indices,data['ice velocity magnitude'][i][j],'->',magnitude
      data['ice velocity azimuth grid'][i][j] = azimuth
      data['ice velocity magnitude'][i][j] = magnitude
    counter[kbc_mask[i,j]] += 1
    kbc_mask[i,j] += 2
  inputfile.close()
  print 'inlets.dat contains',counter[0],'points that are not in kbc.dat'
  print 'inlets.dat contains',counter[1],'points that are in kbc.dat'

########## PART II: CREATE A NETCDF FILE CONTAINING THE RAW DATA ##########

  filename = os.path.join('output','raw.nc')
  print '\nWriting', filename
  if netCDF_module == 'netCDF4':
    netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
  else:
    netCDFfile = NetCDFFile(filename,'w')
  
  netCDFfile.createDimension('x',data['rows columns number of sub parameters'][0][1])
  netCDFfile.createDimension('y', data['rows columns number of sub parameters'][0][0])
  netCDFfile.createVariable('x','f',('x',))[:] = [x[0] for x in data['columns position'][:-1]]
  netCDFfile.createVariable('y','f',('y',))[:] = [x[0] for x in data['rows position'][:-1]]
  netCDFfile.createVariable('mask1',    'i',('y','x'))[:] = data['existency table:']
  netCDFfile.createVariable('azimuth',  'f',('y','x'))[:] = data['ice velocity azimuth grid']
  netCDFfile.createVariable('velocity', 'f',('y','x'))[:] = data['ice velocity magnitude']
  netCDFfile.createVariable('thickness','f',('y','x'))[:] = data['thickness']
  netCDFfile.createVariable('mask2',    'i',('y','x'))[:] = data['reliable velocity obs']
  netCDFfile.createVariable('seabed',   'f',('y','x'))[:] = data['seabed depth']
  netCDFfile.createVariable('mask3',    'i',('y','x'))[:] = data['fake ice shelf region']
  netCDFfile.createVariable('accumulation','f',('y','x'))[:] = data['surface accumulation']
  netCDFfile.createVariable('bbar',        'f',('y','x'))[:] = data['flowlaw']
  netCDFfile.createVariable('temperature', 'f',('y','x'))[:] = data['surface temperature']
  netCDFfile.createVariable('kbc',         'i',('y','x'))[:] = kbc_mask

  netCDFfile.close()
  del(netCDFfile) # remove this variable from the name-space (pycdf might fail if we don't)

########## PART III: CREATE THE NETCDF FILE NEEDED BY GLIMMER ##########

# Read the configuration file to get the number of levels to be used
# Also check that ewn, nsn, dew, and dns are correct
  print '\nReading ross.config'
  configParser = ConfigParser.SafeConfigParser()
  configParser.read('ross.config')
  nx = int(configParser.get('grid','ewn'))
  ny = int(configParser.get('grid','nsn'))
  nz = int(configParser.get('grid','upn'))
  dx = int(configParser.get('grid','dew'))
  dy = int(configParser.get('grid','dns'))

  if nx != 151:
    print 'WARNING: ewn should be set to 151 in ross.config'
  if ny != 115:
    print 'WARNING: nsn should be set to 115 in ross.config'
  if dx != 6822 or dy !=6822:
    print 'WARNING: dew and dns should be set to 6822 in ross.config'

# Put the data into numpy arrays with two extra rows and columns all around
# This reproduces a previous script's output (makerossnc.py)
  mask1        = addBorder(data['existency table:'],dtype='i')
  azimuth      = addBorder(data['ice velocity azimuth grid'],dtype='f')
  velocity     = addBorder(data['ice velocity magnitude'],dtype='f')
  thickness    = addBorder(data['thickness'],dtype='f')
#  mask2        = addBorder(data['reliable velocity obs'],dtype='i')
  seabed       = addBorder(data['seabed depth'],dtype='f',value=5000.0)
  mask3        = addBorder(data['fake ice shelf region'],dtype='i')
#  accumulation = addBorder(data['surface accumulation'],dtype='f')
#  bbar         = addBorder(data['flowlaw'],dtype='f')
#  temperature  = addBorder(data['surface temperature'],dtype='f')
  kbc          = addBorder(kbc_mask,dtype='i')

# Remove any parts of the "fake shelf" that are not in the "existency table"
  if verbose:
    print 'Removing',numpy.sum(numpy.logical_and(mask1==0,mask3!=0)),'points from mask3'
  mask3[mask1==0] = 0

  if mask_llc:
#   Modify mask3 (the "fake shelf") to cut off the lower left corner of the domain
    for ii in range(nx):
      if mask3[2,ii] == 1: break
    if fake_shelf:
      mask3[:,:] = 0
      jj = 0
    else:    
      for jj in range(2,ny):
        if mask3[jj,ii] == 0: break
      jj = jj-2
    for i in range(1,ii):
      for j in range(1,ii+jj-i):
        mask3[j,i] = 1
    for i in range(ii,nx):
      mask3[1,i] = mask1[2,i]
  elif fake_shelf:
    mask3[:,:] = 0

  thickness[mask3==1] = 0

# Set the velocity to zero except where needed as a kinematic boundary condition.
  velocity[kbc == 0] = 0
# Get the components of the velocity vector
  azimuth *= numpy.pi/180
  velocity1 = velocity * numpy.cos(azimuth)
  velocity2 = velocity * numpy.sin(azimuth)

# Create the netCDF input file needed by Glimmer
  filename = os.path.join('output','ross.nc')
  print 'Writing', filename
  if netCDF_module == 'netCDF4':
    netCDFfile = NetCDFFile(filename,'w',format='NETCDF3_CLASSIC')
  else:
    netCDFfile = NetCDFFile(filename,'w')
  netCDFfile.createDimension('time',1)
  netCDFfile.createDimension('x1',nx)
  netCDFfile.createDimension('y1',ny)
  netCDFfile.createDimension('level',nz)
  netCDFfile.createDimension('x0',nx-1) # staggered grid 
  netCDFfile.createDimension('y0',ny-1)
  time  = netCDFfile.createVariable('time','f',('time',))
  x1    = netCDFfile.createVariable('x1','f',('x1',))
  y1    = netCDFfile.createVariable('y1','f',('y1',))
  x0    = netCDFfile.createVariable('x0','f',('x0',))
  y0    = netCDFfile.createVariable('y0','f',('y0',))
  thk   = netCDFfile.createVariable('thk' ,'f',('time','y1','x1'))
  topg  = netCDFfile.createVariable('topg','f',('time','y1','x1'))
  beta  = netCDFfile.createVariable('beta','f',('time','y0','x0'))
  uvelbc = netCDFfile.createVariable('uvelbc','f',('time','level','y0','x0'))
  vvelbc = netCDFfile.createVariable('vvelbc','f',('time','level','y0','x0'))
  uvel = netCDFfile.createVariable('uvel','f',('time','level','y0','x0'))
  vvel = netCDFfile.createVariable('vvel','f',('time','level','y0','x0'))
  kinbcmask = netCDFfile.createVariable('kinbcmask','i',('time','y1','x1'))
  
  time[0] = 0
  x = dx*numpy.arange(nx,dtype='float32')
  y = dx*numpy.arange(ny,dtype='float32')
  x1[:] = x
  y1[:] = y
  x0[:] = dx/2 + x[:-1] # staggered grid
  y0[:] = dy/2 + y[:-1]
  thk [:] = thickness
  topg[:] = -seabed
  beta[:] = numpy.zeros((ny-1,nx-1),dtype='float32')
  uvel[:] = numpy.array(nz*[velocity2[:-1,:-1]])
  vvel[:] = numpy.array(nz*[velocity1[:-1,:-1]])
  mask = numpy.logical_and(velocity==0,numpy.logical_or(mask1==1,mask3==1))
  mask[0,:-1] = True
  mask[:-1,0] = True
  velocity2[mask] = float('NaN')
  velocity1[mask] = float('NaN')
  uvelbc[:] = numpy.array(nz*[velocity2[:-1,:-1]])
  vvelbc[:] = numpy.array(nz*[velocity1[:-1,:-1]])
  kinbcmask[:] = numpy.int32(numpy.where(mask, 0, 1))

  netCDFfile.close()

########## PART IV: RUN GLIMMER ##########

if run_glimmer:
  print '\nRunning',executable,'for the Ross Ice Shelf experiment'
  os.system('echo ross.config'+' | '+executable)

# Clean up by moving extra files written by Glimmer to the 'scratch' subdirectory
print '\nMoving files to the scratch subdirectory: (*) indicates an old file was deleted'
# Look for files with extension 'txt', 'log', or 'nc'
for files in glob.glob('*.txt')+glob.glob('*.log')+glob.glob('*.nc'):
  star = ''
# Delete any files already in scratch with these filenames 
  if files in os.listdir('scratch'):
    os.remove(os.path.join('scratch',files))
    star = '*'
# Move the new files to scratch
  shutil.move(files,'scratch')
  print files+star,
print
