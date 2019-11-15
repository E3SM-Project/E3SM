import netCDF4
import numpy as np
from scipy import spatial
import matplotlib.pyplot as plt
import argparse
plt.switch_backend('agg')

######################################################################################
######################################################################################

def lonlat2xyz(lon,lat):

  R = 6378206.4
  x = R*np.multiply(np.cos(lon),np.cos(lat))
  y = R*np.multiply(np.sin(lon),np.cos(lat))
  z = R*np.sin(lat)

  return x,y,z

######################################################################################
######################################################################################

def create_pointstats_file(mesh_file,stations_files):

  # Read in station locations
  lon = []
  lat = []
  for stations_file in stations_files:
    f = open(stations_file,'r')
    lines = f.read().splitlines()
    for line in lines:
      lon.append(line.split()[0])
      lat.append(line.split()[1])
  
  # Convert station locations
  lon = np.radians(np.array(lon,dtype=np.float32))
  lon_idx, = np.where(lon < 0.0)
  lon[lon_idx] = lon[lon_idx] + 2.0*np.pi 
  lat = np.radians(np.array(lat,dtype=np.float32))
  stations = np.vstack((lon,lat)).T
  #x,y,z = lonlat2xyz(lon,lat)
  #stations = np.vstack((x,y,z)).T
  
  # Read in cell center coordinates
  mesh_nc = netCDF4.Dataset(mesh_file,'r')
  lonCell = np.array(mesh_nc.variables["lonCell"][:])
  latCell = np.array(mesh_nc.variables["latCell"][:])
  meshCells = np.vstack((lonCell,latCell)).T
  #x,y,z = lonlat2xyz(lonCell,latCell)
  #meshCells = np.vstack((x,y,z)).T
  
  # Find nearest cell center to each station
  tree = spatial.KDTree(meshCells)
  d,idx = tree.query(stations)
  
  # Plot the station locations and nearest cell centers
  plt.figure()
  plt.plot(lonCell[idx],latCell[idx],'.')
  plt.plot(lon,lat,'.')
  plt.savefig('station_locations.png')
  
  # Open netCDF file for writing
  data_nc = netCDF4.Dataset('points.nc','w', format='NETCDF3_64BIT_OFFSET')
  
  # Find dimesions
  npts = idx.shape[0]
  ncells = lonCell.shape[0]
  
  # Declare dimensions
  data_nc.createDimension('nCells',ncells)
  data_nc.createDimension('StrLen',64)
  data_nc.createDimension('nPoints',npts)
  
  # Declear variables
  npts = data_nc.dimensions['nPoints'].name
  pnt_ids = data_nc.createVariable('pointCellGlobalID',np.int32,(npts,))
  
  # Set variables
  pnt_ids[:] = idx[:]
  data_nc.close()

######################################################################################
######################################################################################

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser.add_argument('--mesh_file',     help='file that contains the MPAS mesh information')
  parser.add_argument('--station_files', action='append', help='list of files that contain station information')
  args = parser.parse_args()
  
  #mesh_file = 'culled_mesh.nc'
  #stations_files= ['USGS_stations/stations_all.txt','NOAA-COOPS_stations/stations.txt']

  create_pointstats_file(args.mesh_file,args.station_files)
  
