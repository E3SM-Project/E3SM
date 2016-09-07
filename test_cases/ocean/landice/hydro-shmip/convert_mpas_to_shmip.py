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
parser.add_option("-o", "--out", dest="outfile", type='string', help="name of SHMIP format file to create", metavar="FILE")
parser.add_option("-t", "--title", dest="title", type='string', help="string to use for title for this test, e.g. 'hoffman_mpas_A1'", metavar="TITLE")
parser.add_option("-v", "--version", dest="hash", type='string', help="version of MPAS used", metavar="HASH")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'output.nc'
   print 'No file specified.  Attempting to use output.nc'
if not options.icfilename:
   options.icfilename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use '+options.icfilename
if not options.outfile:
   options.outfile = 'SHMIP.nc'
   print 'No output file specified.  Attempting to use '+options.outfile

# Open the file, get needed dimensions
infile = netCDF4.Dataset(options.filename,'r')
icfile = netCDF4.Dataset(options.icfilename,'r')


# Create new file

outfile = netCDF4.Dataset(options.outfile, 'w')

# ============================================
# Create all of the netcdf global attributes
# ============================================
# Do this first as doing it last is slow for big files since adding
# attributes forces the contents to get reorganized.
setattr(outfile, 'title', options.title)
setattr(outfile, 'meshtype', 'unstructured Voronoi')
setattr(outfile, 'dimension', '2D')
setattr(outfile, 'channels_on_edges', 'no')
setattr(outfile, 'institution', 'Matthew Hoffman, Los Alamos National Laboratory')
setattr(outfile, 'source', 'MPAS: '+options.hash)


# ============================================
# Create dimensions
# ============================================
outfile.createDimension('time', None)  # None=unlimited
outfile.createDimension('cellnr', len(infile.dimensions['nCells']))
outfile.createDimension('edgenr', len(infile.dimensions['nEdges']))
outfile.createDimension('dim', 2) # dimension of model domain manifold
outfile.createDimension('layernr', 1) # number of distributed layers
outfile.createDimension('chlayernr', 1) # number of channel layers
outfile.createDimension('maxedges',  len(infile.dimensions['maxEdges']))
outfile.createDimension('nodesperedge',  2)




# ============================================
# Create variables
# ============================================

# DIMENSION VARIABLES
thevar = outfile.createVariable('cellnr', 'i', ('cellnr',))
thevar[:] = infile.variables['indexToCellID'][:]
setattr(thevar, 'description', 'List of global cell IDs.  1-based.')

thevar = outfile.createVariable('edgenr', 'i', ('edgenr',))
thevar[:] = infile.variables['indexToEdgeID'][:]
setattr(thevar, 'description', 'List of global edge IDs.  1-based.')

thevar = outfile.createVariable('dim', 'i', ('dim',))
thevar[:] = (1,2)

thevar = outfile.createVariable('time', 'd', ('time',))
thevar[0] = infile.variables['daysSinceStart'][-1] * 3600.0 * 24.0
setattr(thevar, 'units', 's')

thevar = outfile.createVariable('layernr', 'i', ('layernr',))
thevar[:] = (1,)
setattr(thevar, 'layers_description', 'macroporous')

thevar = outfile.createVariable('chlayernr', 'i', ('chlayernr',))
thevar[:] = (1,)
setattr(thevar, 'channel_layers_description', 'cross sectional area')

thevar = outfile.createVariable('xy', 'd', ('dim','cellnr'))
thevar[0,:] = infile.variables['xCell'][:]
thevar[1,:] = infile.variables['yCell'][:]
setattr(thevar, 'units', 'm')
setattr(thevar, 'description', 'coordinates of cell centers')


# CONNECTIVITY
thevar = outfile.createVariable('cellconnect', 'i', ('maxedges', 'cellnr',))
thevar[:] = np.transpose(infile.variables['cellsOnCell'][:])
setattr(thevar, 'description', 'List of cells that neighbor each cell.  Value of 0 indicates this cell does not have that many edges.')

thevar = outfile.createVariable('edgeconnect', 'i', ('nodesperedge', 'edgenr',))
thevar[:] = np.transpose(infile.variables['cellsOnEdge'][:])
setattr(thevar, 'description', 'List of cells that neighbor each edge.')


# HYDRO VARIABLES
thevar = outfile.createVariable('H', 'd', ('cellnr',))
thevar[:] = icfile.variables['thickness'][-1,0]
setattr(thevar, 'description', 'ice thickness')
setattr(thevar, 'units', 'm')

thevar = outfile.createVariable('B', 'd', ('cellnr',))
thevar[:] = icfile.variables['bedTopography'][-1,0]
setattr(thevar, 'description', 'bed elevation')
setattr(thevar, 'units', 'm')

thevar = outfile.createVariable('N', 'd', ('time', 'cellnr',))
thevar[0,:] = infile.variables['effectivePressure'][-1,0]
setattr(thevar, 'description', 'effective pressure')
setattr(thevar, 'units', 'Pa')

thevar = outfile.createVariable('hs', 'd', ('time', 'cellnr',))
thevar[0,:] = infile.variables['waterThickness'][-1,0]
setattr(thevar, 'description', 'effective layer thickness')
setattr(thevar, 'units', 'm')

thevar = outfile.createVariable('Cs', 'd', ('time', 'edgenr',))
thevar[0,:] = infile.variables['channelArea'][-1,0]
setattr(thevar, 'description', 'channel cross-sectional area')
setattr(thevar, 'units', 'm^2')

thevar = outfile.createVariable('Qs', 'd', ('time', 'edgenr',))
thevar[0,:] = np.absolute(infile.variables['channelDischarge'][-1,0])
setattr(thevar, 'description', 'channel discharge (absolute value)')
setattr(thevar, 'units', 'm^3/s')

thevar = outfile.createVariable('qs', 'd', ('time', 'edgenr',))
thevar[0,:] = infile.variables['waterFlux'][-1,0]
setattr(thevar, 'description', 'discharge in the distributed layers (absolute value)')
setattr(thevar, 'units', 'm^2/s')


outfile.close()

