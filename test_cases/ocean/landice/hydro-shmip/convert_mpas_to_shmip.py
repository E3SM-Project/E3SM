#!/usr/bin/env python
'''
Convert MPAS hydro output to the SHMIP format regquired by SHMIP
Details here: http://shmip.bitbucket.org/technical-instructions.html
'''
import netCDF4
import numpy as np
import sys

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="MPAS output file", metavar="FILE")
parser.add_option("-i", "--icfile", dest="icfilename", type='string', help="MPAS initial condition file", metavar="FILE")
parser.add_option("-s", "--scenario", dest="scenario", type='string', help="name of SHMIP scenario this run corresponds to, e.g., A1", metavar="SCENARIO")
parser.add_option("-t", "--title", dest="title", type='string', help="string to use for title for this test, e.g. 'hoffman_mpas_A1'", metavar="TITLE")
parser.add_option("-v", "--version", dest="hash", type='string', help="version of MPAS used", metavar="HASH")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'output.nc'
   print 'No file specified.  Attempting to use output.nc'
if not options.icfilename:
   options.icfilename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use '+options.icfilename
if not options.scenario:
   print 'ERROR: Scenario has not been specified.  Please include this with the -s option.'
   sys.exit()
if not options.title:
   options.title = 'hoffman_mpas_'+options.scenario
   print 'No title specified.  Using: '+options.title
if not options.hash:
   options.hash = 'NA'
   print 'Not version specified.  Using: '+options.hash

# Open the file, get needed dimensions
infile = netCDF4.Dataset(options.filename,'r')
icfile = netCDF4.Dataset(options.icfilename,'r')


# Create new file

outfilename = "{}_mhof.nc".format(options.scenario)
outfile = netCDF4.Dataset(outfilename, 'w')

# ============================================
# Create all of the netcdf global attributes
# ============================================
# Do this first as doing it last is slow for big files since adding
# attributes forces the contents to get reorganized.
setattr(outfile, 'title', options.title)
setattr(outfile, 'meshtype', 'unstructured Voronoi')
setattr(outfile, 'dimension', '2D')
setattr(outfile, 'channels_on_edges', 'yes')
setattr(outfile, 'institution', 'Matthew Hoffman, Los Alamos National Laboratory')
setattr(outfile, 'source', 'MPAS: '+options.hash)
setattr(outfile, 'references', 'http://shmip.bitbucket.io/')


# ============================================
# Create dimensions
# ============================================
outfile.createDimension('time', None)  # None=unlimited
outfile.createDimension('dim', 2) # spatial dimensions
outfile.createDimension('n_nodes_ch', 2) # how many nodes to make a channel edge.  Fixed(?) at 2
outfile.createDimension('index1', len(infile.dimensions['nCells']))
outfile.createDimension('index2', len(infile.dimensions['nEdges']))
outfile.createDimension('index_ch', len(infile.dimensions['nEdges']))


# ============================================
# Account for time depending on test
# ============================================
ntIn = len(infile.dimensions['Time'])
if options.scenario[0] in ('C','D'):
   times = np.arange(ntIn)  # use list of all time indices.  (Assuming file ONLY contains the required time levels (i.e. annual file with daily output for test D; daily file with hourly output for test C)
else:
   times = np.array([ntIn,])  # only index of final time
print "Using time indices:", times



# ============================================
# Create variables
# ============================================

# DIMENSION VARIABLES

thevar = outfile.createVariable('time', 'd', ('time',))
thevar[times-times.min()] = (infile.variables['daysSinceStart'][times] - infile.variables['daysSinceStart'][times[0]])* 3600.0 * 24.0
setattr(thevar, 'units', 's')
setattr(thevar, 'long_name', 'time')

thevar = outfile.createVariable('coords1', 'd', ('dim','index1'))
thevar[0,:] = infile.variables['xCell'][:]
thevar[1,:] = infile.variables['yCell'][:]
setattr(thevar, 'units', 'm')
setattr(thevar, 'long_name', 'node coordinates')

thevar = outfile.createVariable('coords2', 'd', ('dim','index2'))
thevar[0,:] = infile.variables['xEdge'][:]
thevar[1,:] = infile.variables['yEdge'][:]
setattr(thevar, 'units', 'm')
setattr(thevar, 'long_name', 'cell midpoint coordinates')

thevar = outfile.createVariable('coords_ch', 'd', ('dim','index_ch'))
thevar[0,:] = infile.variables['xEdge'][:]
thevar[1,:] = infile.variables['yEdge'][:]
setattr(thevar, 'units', 'm')
setattr(thevar, 'long_name', 'channel midpoint coordinates')

thevar = outfile.createVariable('connect_ch', 'i', ('n_nodes_ch', 'index_ch',))
thevar[:] = np.transpose(infile.variables['cellsOnEdge'][:]-1)  # convert from 1-based to 0-based
setattr(thevar, 'units', '')
setattr(thevar, 'long_name', 'channel connectivity')



# GEOMETRY/HYDRO VARIABLES

thevar = outfile.createVariable('H', 'd', ('index1',))
thevar[:] = icfile.variables['thickness'][-1,:]
setattr(thevar, 'long_name', 'ice thickness')
setattr(thevar, 'units', 'm')

thevar = outfile.createVariable('B', 'd', ('index1',))
thevar[:] = icfile.variables['bedTopography'][-1,:]
setattr(thevar, 'long_name', 'bed elevation')
setattr(thevar, 'units', 'm')

thevar = outfile.createVariable('N', 'd', ('time', 'index1',))
thevar[times-times.min(),:] = infile.variables['effectivePressure'][times,:]
setattr(thevar, 'long_name', 'effective pressure')
setattr(thevar, 'units', 'Pa')

# Optional additional variables

thevar = outfile.createVariable('h', 'd', ('time', 'index1',))
thevar[times-times.min(),:] = infile.variables['waterThickness'][times,:]
setattr(thevar, 'long_name', 'water sheet thickness')
setattr(thevar, 'units', 'm')

# not used: stored water effective layer thickness

# Edge variables

thevar = outfile.createVariable('q', 'd', ('time', 'index2',))
thevar[times-times.min(),:] = np.absolute(infile.variables['waterFlux'][times,:])
setattr(thevar, 'long_name', 'water sheet discharge')
setattr(thevar, 'units', 'm^2/s')

# Channel variables

thevar = outfile.createVariable('S', 'd', ('time', 'index_ch',))
thevar[times-times.min(),:] = infile.variables['channelArea'][times,:]
setattr(thevar, 'long_name', 'channel cross-sectional area')
setattr(thevar, 'units', 'm^2')

thevar = outfile.createVariable('Q', 'd', ('time', 'index_ch',))
thevar[times-times.min(),:] = np.absolute(infile.variables['channelDischarge'][times,:])
setattr(thevar, 'long_name', 'channel discharge')
setattr(thevar, 'units', 'm^3/s')


outfile.close()

print '\nConversion complete.  Written to: '+outfilename
