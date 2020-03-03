import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import glob
import pprint
import datetime
import os
import yaml
import subprocess
import argparse
import write_forcing_file
plt.switch_backend('agg')

##################################################################################################
##################################################################################################

if __name__ == '__main__':

  # This creates a "dummy" time varying forcing file
  # with zero wind zero atmospheric pressure perturbation
  # for the tidal spinup run.
  #
  # The tidal spinup is run using this "dummy" atmospheric forcing
  # because the time varying atmospheric forcing for the 
  # forward run requires information in the restart file.
  # The inclusion of this additional information in the restart
  # file is trigged by the use of time varying atmospheric forcing
  # in the tidal spinup.

  parser = argparse.ArgumentParser()
  parser.add_argument('--start_time')
  parser.add_argument('--spinup_length')
  args = parser.parse_args()

  # Files to interpolate to/from
  grid_file = './mesh.nc'
  forcing_file = 'spinup_atmospheric_forcing.nc'

  # Setup timestamps 
  # (3 time snaps are needed because new data will be read in at the end of the simulation)
  dtformat = '%Y-%m-%d_%H:%M:%S'
  start_time = datetime.datetime.strptime(args.start_time,dtformat)
  spinup_length = float(args.spinup_length)
  xtime = []
  xtime.append(args.start_time+45*' ')
  next_time = start_time + datetime.timedelta(days=spinup_length)
  xtime.append(datetime.datetime.strftime(next_time,dtformat)+45*' ')
  next_time = next_time + datetime.timedelta(days=spinup_length)
  xtime.append(datetime.datetime.strftime(next_time,dtformat)+45*' ')
  xtime = np.array(xtime,'S64')
  print(xtime)

  # Get grid from grid file
  grid_nc = netCDF4.Dataset(grid_file,'r')
  lon_grid = grid_nc.variables['lonCell'][:]
  ncells = lon_grid.size

  # Initialize atmospheric forcing fields
  u_data = np.zeros((3,ncells))
  v_data = np.zeros((3,ncells))
  p_data = np.zeros((3,ncells)) + 101325.0
  print(p_data.shape)

  # Write to NetCDF file
  subprocess.call(['rm',forcing_file])
  write_forcing_file.write_to_file(forcing_file,u_data,'windSpeedU',xtime)
  write_forcing_file.write_to_file(forcing_file,v_data,'windSpeedV',xtime)
  write_forcing_file.write_to_file(forcing_file,p_data,'atmosPressure',xtime)

