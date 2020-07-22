# Author: Zhendong Cao
# Date: July 2020
# Description: This python script is used to generate vegetationInfo_Delaware.nc
#              by reading and interpolating the data of vegetation properties and
#              Manning's roughness coefficients from tif files to the MPAS-O grid file.

####################################################################################################
#################################### function to read tif files ####################################
def read_tif(filename):

    gtif = gdal.Open(filename)
    if gtif is None:
        print ('Unable to open tif file')
        sys.exit(1)
    band1 = gtif.GetRasterBand(1)
    geotransform = gtif.GetGeoTransform()
    data = band1.ReadAsArray()

    # set negative vegetation data to zero
    data[data<0.0]=0.0

    # convert gettransform info to lat/lon
    orig_lon, lon_res, a, orig_lat, b, lat_res = geotransform
    nlat, nlon = data.shape
    lat=[]
    lon=[]

    for i in range(nlat):
        lat_i = orig_lat + i * lat_res
        lat.append(lat_i)

    for j in range(nlon):
        lon_j = orig_lon + j * lon_res
        lon.append(lon_j)

    return data, lat, lon

def mesh_tri(grid_file):

  ncmesh = Dataset(grid_file,'r')
  lon_mesh = np.rad2deg(ncmesh.variables['lonCell'][:])-360.0
  #lon_mesh = np.rad2deg(np.mod(ncmesh.variables['lonCell'][:] + np.pi, 2.0*np.pi) - np.pi)
  lat_mesh = np.rad2deg(ncmesh.variables['latCell'][:])
  nEdgesOnCell = ncmesh.variables['nEdgesOnCell'][:]
  cellsOnCell = ncmesh.variables['cellsOnCell'][:,:]
 
  # Triangulate cells 
  triangles = Triangulation(lon_mesh,lat_mesh)

  # Compute triangulation mask (needs to be vectorized)
  mask = np.array(np.zeros((triangles.triangles.shape[0],)),dtype=bool)
  ntri = triangles.neighbors.shape[0]
  for i in range(ntri):
    k = 0
    for j in range(3):
      n = triangles.triangles[i,j]
      if nEdgesOnCell[n] != np.where(cellsOnCell[n,:] != 0)[0].size:  # Mask triangles
        k = k +1                                                      # containing
    if k == 3:                                                        # 3 boundary
      mask[i] = True                                                  # cells
  triangles.set_mask(mask)


  return triangles

######################## function to interpolate tif data to MPAS-O grid ###########################
def interpolate_data_to_grid(data_file,grid_file,fill=0):

    # get tif data info
    data, lat_data, lon_data = read_tif(data_file)

    # get MPAS-O grid data info
    grid_nc = Dataset(grid_file, 'r')
    lon_grid = grid_nc.variables['lonCell'][:]*180.0/np.pi-360
    lat_grid = grid_nc.variables['latCell'][:]*180.0/np.pi
    grid_points = np.column_stack((lon_grid, lat_grid))
    ncells = lon_grid.size
    interp_data = np.zeros((ncells))

    # interpolate tif data onto new grid
    interpolator = grid_interpolator((lon_data,lat_data[::-1]),np.flipud(data).T,
                                     method='nearest',
                                     bounds_error=False,fill_value=fill)
    interp_data = interpolator(grid_points)

    return lon_data,lat_data,data,lon_grid,lat_grid,interp_data

############################### function to plot interpolated data #################################
def plot_interp_data(lon_data,lat_data,data,lon_grid,lat_grid,interp_data,tri):

    # plot original data
    fig, ax = plt.subplots(1,2, figsize=(12,5))
    levels = np.linspace(np.amin(data),np.amax(data),10)
    lonn_data, latt_data = np.meshgrid(lon_data, lat_data)
    cf = ax[0].contourf(lonn_data, latt_data, data, levels=levels)
    ax[0].set_title('original data')
    cbar = fig.colorbar(cf,ax=ax[0])

    # plot interpolated data
    lon_max = np.max(lon_data) + 5.0
    lon_min = np.min(lon_data) - 5.0
    lat_max = np.max(lat_data) + 5.0
    lat_min = np.min(lat_data) - 5.0
    levels = np.linspace(np.amin(interp_data),np.amax(interp_data),10)
    cf = ax[1].tricontourf(tri, interp_data, levels=levels)
    ax[1].set_xlim([lon_min,lon_max])
    ax[1].set_ylim([lat_min,lat_max])
    ax[1].set_title('interpolated data')
    cbar = fig.colorbar(cf,ax=ax[1])

####################################################################################################
#########################################  MAIN CODE  ##############################################
if __name__ == '__main__':
  import yaml
  import gdal,osr
  import os
  import sys
  import pprint
  import numpy as np
  import matplotlib.pyplot as plt
  from matplotlib.tri import Triangulation
  from scipy.interpolate import RegularGridInterpolator as grid_interpolator
  from netCDF4 import Dataset

################################ read config file  #################################################
  pwd = os.getcwd()
  f = open(pwd+'/vegetation_tif_to_mpas.config')
  cfg = yaml.load(f, Loader=yaml.Loader)
  pprint.pprint(cfg)

  vegMask_file = cfg['input_location'] + cfg['vegetation_mask_file']
  vegHeight_file = cfg['input_location'] + cfg['vegetation_height_file']
  vegDensity_file = cfg['input_location'] + cfg['vegetation_density_file']
  vegDiameter_file = cfg['input_location'] + cfg['vegetation_diameter_file']
  Manning_file = cfg['input_location'] + cfg['Delaware_Manning_file']
  grid_file = cfg['mpas_grid_file']
  bottomDrag_file = cfg['bottomDrag_file']
  figplot = cfg['figplot']
  output_file = cfg['output_file']

  tri = mesh_tri(grid_file)
############################### data interpolation and plot ########################################
  ##vegetation Mask
  lon_veg, lat_veg, vegMask, lon_grid, lat_grid, interp_vegMask = interpolate_data_to_grid(vegMask_file, grid_file)
  vegMask[vegMask==0]=0
  vegMask[vegMask>0]=1
  interp_vegMask[interp_vegMask==0]=0
  interp_vegMask[interp_vegMask>0]=1
  interp_vegMask = interp_vegMask.astype(int)
  if (figplot):
    plot_interp_data(lon_veg, lat_veg, vegMask, lon_grid, lat_grid, interp_vegMask, tri)
    plt.savefig('vegetation_mask.png')

  ##vegHeight
  lon_veg, lat_veg, vegHght, lon_grid, lat_grid, interp_vegHght = interpolate_data_to_grid(vegHeight_file, grid_file)
  if (figplot):
    plot_interp_data(lon_veg, lat_veg, vegHght, lon_grid, lat_grid, interp_vegHght, tri)
    plt.savefig('vegetation_height.png')

  ##vegDensity
  lon_veg, lat_veg, vegDens, lon_grid, lat_grid, interp_vegDens = interpolate_data_to_grid(vegDensity_file, grid_file)
  if (figplot):
    plot_interp_data(lon_veg, lat_veg, vegDens, lon_grid, lat_grid, interp_vegDens, tri)
    plt.savefig('vegetation_density.png')

  ##vegDiameter
  lon_veg, lat_veg, vegDiam, lon_grid, lat_grid, interp_vegDiam = interpolate_data_to_grid(vegDiameter_file, grid_file)
  if (figplot):
    plot_interp_data(lon_veg, lat_veg, vegDiam, lon_grid, lat_grid, interp_vegDiam, tri)
    plt.savefig('vegetation_diameter.png')

  ##Manning
  forcing_nc = Dataset(bottomDrag_file,'r')
  manning_init = forcing_nc.variables['bottomDrag'][:]
  lon_manning, lat_manning, manning, lon_grid, lat_grid, interp_manning = interpolate_data_to_grid(Manning_file, grid_file, fill=-9999)
  interp_manning[interp_manning < 1e-5] = manning_init[interp_manning < 1e-5]
  if (figplot):
    plot_interp_data(lon_manning, lat_manning, manning, lon_grid, lat_grid, manning_init, tri)
    plt.savefig('manning_init.png')
    plot_interp_data(lon_manning, lat_manning, manning, lon_grid, lat_grid, interp_manning, tri)
    plt.savefig('vegetation_manning.png')

################################# output file generation ###########################################
  if os.path.exists(output_file):
    os.remove(output_file)

  newfile = Dataset(output_file,'w', format='NETCDF3_64BIT_OFFSET')

  # define dimension
  newfile.createDimension('nCells',len(lat_grid))

  # create variables
  veg_mask = newfile.createVariable('vegetationMask', np.int32, ('nCells'))
  veg_mask.standard_name='vegetation_mask'
  veg_mask.units =''

  veg_diam = newfile.createVariable('vegetationDiameter', np.float64, ('nCells'))
  veg_diam.standard_name='stem diameter of each vegetation shoot'
  veg_diam.units ='m'

  veg_hght = newfile.createVariable('vegetationHeight', np.float64, ('nCells'))
  veg_hght.standard_name='stem height of each vegetation shoot'
  veg_hght.units ='m'

  veg_dens = newfile.createVariable('vegetationDensity', np.float64, ('nCells'))
  veg_dens.standard_name='stem numbers per unit area'
  veg_dens.units ='m^{-2}'

  Delaware_manning = newfile.createVariable('bottomDrag', np.float64, ('nCells'))
  Delaware_manning.standard_name='bottom Manning roughness coefficient'
  Delaware_manning.units ='s/m^{1/3}'

  # assign variable values
  veg_mask[:] = interp_vegMask[:]
  veg_diam[:] = interp_vegDiam[:]
  veg_hght[:] = interp_vegHght[:]
  veg_dens[:] = interp_vegDens[:]
  Delaware_manning[:] = interp_manning[:]

  # close file
  newfile.close()

##!---------- end of script ---------
