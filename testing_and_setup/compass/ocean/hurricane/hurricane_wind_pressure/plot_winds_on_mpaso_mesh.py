# Author: Steven Brus
# Date: April, 2020
# Description: Plots syntetic wind/pressure timeseries on MPAS-O mesh

import netCDF4
import matplotlib.pyplot as plt
import numpy as np
import os
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
plt.switch_backend('agg')
cartopy.config['pre_existing_data_dir'] = \
    os.getenv('CARTOPY_DIR', cartopy.config.get('pre_existing_data_dir'))

#######################################################################
#######################################################################

def plot_data(lon_grid,lat_grid,data,var_label,var_abrev,time):

  fig = plt.figure()
  ax1 = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
  levels = np.linspace(np.amin(data),np.amax(data),100)
  cf = ax1.tricontourf(lon_grid,lat_grid,data,levels=levels,transform=ccrs.PlateCarree())
  ax1.set_extent([0, 359.9, -90, 90], crs=ccrs.PlateCarree())
  ax1.add_feature(cfeature.LAND, zorder=100)
  ax1.add_feature(cfeature.LAKES, alpha=0.5, zorder=101)
  ax1.add_feature(cfeature.COASTLINE, zorder=101)
  ax1.set_title('interpolated data '+time.strip())
  cbar = fig.colorbar(cf,ax=ax1)
  cbar.set_label(var_label)
  
  # Save figure
  fig.tight_layout()
  fig.savefig(var_abrev+'_'+str(i).zfill(4)+'.png',box_inches='tight')
  plt.close()

#######################################################################
#######################################################################

if __name__ == '__main__':
  
  grid_file = 'mesh.nc'
  data_file = 'out.nc'

  grid_nc = netCDF4.Dataset(grid_file,'r')
  lon_grid = grid_nc.variables['lonCell'][:]*180.0/np.pi
  lat_grid = grid_nc.variables['latCell'][:]*180.0/np.pi

  data_nc = netCDF4.Dataset(data_file,'r')
  u_data = data_nc.variables['windSpeedU'][:]
  v_data = data_nc.variables['windSpeedV'][:]
  p_data = data_nc.variables['atmosPressure'][:]
  xtime = data_nc.variables['xtime'][:]

  for i in range(u_data.shape[0]-1):

      print('Plotting vel: '+str(i))

      data = np.sqrt(np.square(u_data[i,:]) + np.square(v_data[i,:]))
      time_ls = [x.decode("utf-8") for x in xtime[i]]
      time = ''.join(time_ls)
      plot_data(lon_grid,lat_grid,data,'velocity magnitude','vel',time)
      plot_data(lon_grid,lat_grid,p_data[i,:],'atmospheric pressure','pres',time)
