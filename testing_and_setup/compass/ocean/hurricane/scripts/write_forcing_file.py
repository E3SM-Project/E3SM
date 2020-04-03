# Author: Steven Brus
# Date: April, 2020
# Description: This function writes time-varying forcing data to an input file for the model run.

import os
import numpy as np
import netCDF4

##################################################################################################
##################################################################################################

def write_to_file(filename,data,var,xtime):

  if os.path.isfile(filename):
    data_nc = netCDF4.Dataset(filename,'a', format='NETCDF3_64BIT_OFFSET')
  else:
    data_nc = netCDF4.Dataset(filename,'w', format='NETCDF3_64BIT_OFFSET')

    # Find dimesions
    ncells = data.shape[1]
    nsnaps = data.shape[0]

    # Declare dimensions
    data_nc.createDimension('nCells',ncells)
    data_nc.createDimension('StrLen',64)
    data_nc.createDimension('Time',None)

    # Create time variable
    time = data_nc.createVariable('xtime','S1',('Time','StrLen'))
    time[:,:] = netCDF4.stringtochar(xtime)

  # Set variables
  data_var = data_nc.createVariable(var,np.float64,('Time','nCells'))
  data_var[:,:] = data[:,:]
  data_nc.close()

##################################################################################################
##################################################################################################
