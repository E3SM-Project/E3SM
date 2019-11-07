#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, unicode_literals
import netCDF4
import numpy as np

f = netCDF4.Dataset('SMB_NorESM_2.6_1995-2100.anomaly.nc','r+')

nt = len(f.dimensions['Time'])

# It seems much faster for large files to use NCO to add xtime, like:
# ncap2 -s 'defdim("StrLen",64); xtime[$Time,$StrLen]=" "' TF_NorESM_2.6_1995-2100.nc temp.nc
StrLen=64
if 'xtime' in f.variables:
    xtime = f.variables['xtime']
else:
    if not 'StrLen' in f.dimensions:
       print("adding StrLen")
       f.createDimension('StrLen', StrLen)
    print("adding xtime")
    xtime = f.createVariable('xtime','c', ('Time', 'StrLen'))


startYear = 1995
for i in range(nt):

    time_string = f.variables['xtime'][i,:]

    time_string = "{:04}-01-01_00:00:00".format(i+startYear)

    print(time_string)

    f.variables['xtime'][i,:] = list(time_string.ljust(StrLen))


f.close()

