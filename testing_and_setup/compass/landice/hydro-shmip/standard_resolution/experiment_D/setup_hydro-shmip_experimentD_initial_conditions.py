#!/usr/bin/env python
'''
Generate initial conditions for hydro-SHMIP land ice test case D
Details here: http://shmip.bitbucket.org/
'''

from netCDF4 import Dataset as NetCDFFile
import netCDF4
import numpy as np
import sys
import shutil

# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file to setup", metavar="FILE")
parser.add_option("-a", dest="afile", type='string', help="restart file from test A1 to use to set up this test", metavar="FILE")
#parser.add_option("-n", "--number", dest="number", type='int', help="test variant to set up, 1-5", metavar="NUMBER")  # NOT USED
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'

if not options.afile:
   sys.exit("Error: A restart file from test A1 is required to set up this test.  Specify with -a")

# copy the restart file to be the new input file
shutil.copyfile(options.afile, options.filename)

# Open the file, get needed dimensions
gridfile = NetCDFFile(options.filename,'r+')
StrLen = len(gridfile.dimensions['StrLen'])
gridfile.variables['xtime'][0,:] = netCDF4.stringtoarr('0000-01-01_00:00:00'.ljust(StrLen), StrLen)
gridfile.variables['simulationStartTime'][:] = netCDF4.stringtoarr('0000-01-01_00:00:00'.ljust(StrLen), StrLen)

# modify melt inputs - basalMeltInput remains the same as A1
gridfile.variables['basalMeltInput'][0,:] = 7.93e-11 * 1000.0  # Put background input here
gridfile.close()

print 'Successfully added initial conditions to: ', options.filename

