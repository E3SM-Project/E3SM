import netCDF4 
import numpy as np
import hurricane_model as Hurricane
import structures as Geogrid
import winds_io as WindModel
import matplotlib.pyplot as plt
import datetime

def write_netcdf(filename: str, curr_hurricane: Hurricane, grid: Geogrid, winds: WindModel):
    # http://unidata.github.io/netcdf4-python/#section1
    rootgrp = netCDF4.Dataset(filename, "w", format="NETCDF3_64BIT_OFFSET")

    # Declare dimensions
    rootgrp.createDimension('nCells',grid.ncells)
    rootgrp.createDimension('StrLen',64)
    rootgrp.createDimension('Time',None)

    # Declare variables
    time = rootgrp.dimensions['Time'].name
    ncells = rootgrp.dimensions['nCells'].name
    time_var = rootgrp.createVariable('xtime','S1',('Time','StrLen'))
    u_var = rootgrp.createVariable('windSpeedU',np.float64,(time,ncells))
    v_var = rootgrp.createVariable('windSpeedV',np.float64,(time,ncells))
    pres_var = rootgrp.createVariable('atmosPressure',np.float64,(time,ncells))

    # Format time
    ref_date = curr_hurricane[0].ref_time 
    xtime = []
    for it in range(0,len(curr_hurricane)-1):
      t = curr_hurricane[it].time
      date = ref_date + datetime.timedelta(hours=np.float64(t))
      xtime.append(date.strftime('%Y-%m-%d_%H:%M:%S'+45*' '))
    xtime = np.asarray(xtime)
    xtime_list = []
    for t in xtime:
      xtime_list.append(list(t))
    time_var[:] = xtime_list

    # Assign variables
    kmh_to_mps = 0.277778
    mbar_to_pa = 100.0
    for it in range(0, len(curr_hurricane)-1):
        u_var[it, :] = winds[it].u * kmh_to_mps
        v_var[it, :] = winds[it].v * kmh_to_mps
        pres_var[it, :] = winds[it].pressure_profile * mbar_to_pa

    # Close
    rootgrp.close()
