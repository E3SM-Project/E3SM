#!/usr/bin/env python
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import xarray as xr
import netCDF4
import datetime

forcing = xr.open_dataset('forcing.nc')

# build xtime
time = np.arange(26)
ref_date = datetime.datetime.strptime('1901-01-01_00:00:00','%Y-%m-%d_%H:%M:%S')
xtime = []
for t in time:
  date = ref_date + datetime.timedelta(hours=np.float64(t))
  xdate = date.strftime('%Y-%m-%d_%H:%M:%S')+45*' '
  xdate = '0001'+xdate[4:]                                  # Change year because strftime has issues with years <1900
  xtime.append(xdate)

scaling = np.linspace(0,1.0,len(xtime)).T

# build fields
def make_data_array(values, scaling, forcing):
    #adjusted = np.sign(values)*np.sqrt(scaling[np.newaxis,:].T*1e3*abs(values))
    adjusted = np.repeat(100.0*values,len(xtime),axis=0)*scaling[:,np.newaxis]
    return adjusted

windSpeedU = make_data_array(forcing.windStressZonal.values, scaling, forcing)
windSpeedV = make_data_array(forcing.windStressMeridional.values, scaling, forcing)
atmosphericPressure = windSpeedU*0.0 + 101325.0

ncds = netCDF4.Dataset('atmospheric_forcing.nc', 'w', format='NETCDF3_64BIT_OFFSET')

ncds.createDimension('nCells', len(forcing.nCells))
ncds.createDimension('StrLen', 64)
ncds.createDimension('Time', None)

time = ncds.createVariable('xtime','S1', ('Time', 'StrLen'))
time[:] = netCDF4.stringtochar(np.asarray(xtime, dtype='S64'))

time = ncds.dimensions['Time'].name
ncells = ncds.dimensions['nCells'].name
ncds.createVariable('windSpeedU', np.float64,(time, ncells))
ncds.createVariable('windSpeedV', np.float64,(time, ncells))
ncds.createVariable('atmosPressure', np.float64,(time, ncells))

ncds.variables['windSpeedU'][:,:] = windSpeedU[:,:]
ncds.variables['windSpeedV'][:,:] = windSpeedV[:,:]
ncds.variables['atmosPressure'][:,:] = atmosphericPressure[:,:]

ncds.close()
