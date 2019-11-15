#!/usr/bin/env python
from netCDF4 import Dataset
from optparse import OptionParser

parser = OptionParser()

options, args = parser.parse_args()
fileName=args[0]

nc = Dataset(fileName,'r+')
nc.setncattr('is_periodic', 'NO')
nc.close()

