#!/usr/bin/env python3
#---------------------------------------------------------------------------------------------------
'''
This is a replacement for the legacy gen_domain tool created for CESM. 
Most legacy functionality is reproduced, with the notable exception of 
the pole point latitude adjustment needed for the CESM FV grid.

Created April, 2024 by Walter Hannah (LLNL) 
'''
#---------------------------------------------------------------------------------------------------
'''
The map file used to generate the atmsrf files can be created a few different ways.
For a typical E3SM configuration we recommend using a conservative, monotone map.
Here is an example command that can be used to generate one as of NCO version 5.2.2

  SRC_GRID=${DIN_LOC_ROOT}/../mapping/grids/1x1d.nc
  DST_GRID=${GRID_ROOT}/ne${NE}pg2_scrip.nc
  MAP_FILE=${MAP_ROOT}/map_1x1_to_ne${NE}pg2_traave.nc
  ncremap -a traave --src_grd=${SRC_GRID} --dst_grd=${DST_GRID} --map_file=${MAP_FILE}

'''
#---------------------------------------------------------------------------------------------------
import datetime, os, numpy as np, xarray as xr, numba, itertools, subprocess as sp
user, host = os.getenv('USER'), os.getenv('HOST')
source_code_meta = 'generate_domain_E3SM.py'
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
usage = '''
python mkatmsrffile.py  --map_file <map_file_path> 
                        --dst_grid <atm_grid_name>
                        [--output-root <path>]
                        [--date-stamp <date string>]

Purpose:
  This tool generates a file needed to prescribe atmospheric dry deposition of aerosols at the surface.

Environment
  
  This tool requires a few special packages, such as xarray, numba, and itertools.
  These are all included in the E3SM unified environment:
  https://e3sm.org/resources/tools/other-tools/e3sm-unified-environment/

  Otherwise a simple conda environment can be created:
  conda create --name example_env --channel conda-forge xarray numpy numba scikit-learn netcdf4

The following output file is created:

  atmsrf_<dst_grid>_<date_stamp>.nc
    ????

'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--map_file',
                  dest='map_file',
                  default=None,
                  help='Mapping file from a 1x1 degree grid to the target atmosphere grid')
parser.add_option('--vegetation_file',
                  dest='vegetation_file',
                  default=None,
                  help='Path to file containing area fractions of surface types and plant functional type (PFT) data (regrid_vegetation.nc)')
parser.add_option('--soil_water_file',
                  dest='soil_water_file',
                  default=None,
                  help='Path to file containing soil water data (clim_soilw.nc)')
parser.add_option('--dst_grid',
                  dest='dst_grid',
                  default=None,
                  help='destination atmosphere grid name')
parser.add_option('--output_root',
                  dest='output_root',
                  default='./',
                  help='Output path for domain files')
parser.add_option('--date-stamp',
                  dest='date_stamp',
                  default=None,
                  help='Creation date stamp for domain files')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
def main():
  #-----------------------------------------------------------------------------
  # check for valid input arguments
  if opts.map_file is None: 
    raise ValueError(f'{clr.RED}map_file was not specified{clr.END}')
  if opts.vegetation_file is None: 
    raise ValueError(f'{clr.RED}vegetation_file was not specified{clr.END}')
  if opts.soil_water_file is None: 
    raise ValueError(f'{clr.RED}soil_water_file was not specified{clr.END}')
  if opts.dst_grid is None: 
    raise ValueError(f'{clr.RED}dst_grid was not specified{clr.END}')
  #-----------------------------------------------------------------------------
  # check existence of input files and output path
  if not os.path.exists(opts.vegetation_file) :
     raise ValueError(f'{clr.RED}Vegetaion file does not exist:{clr.END} {opts.vegetation_file}')
  if not os.path.exists(opts.soil_water_file) :
     raise ValueError(f'{clr.RED}Soil water file does not exist:{clr.END} {opts.soil_water_file}')
  if not os.path.exists(opts.output_root) :
    raise ValueError(f'{clr.RED}Output root path does not exist:{clr.END} {opts.output_root}')

  #-----------------------------------------------------------------------------
  # Set date stamp for file name

  if opts.date_stamp is None:
    cdate = datetime.datetime.now().strftime('%Y%m%d')
  else:
    cdate = opts.date_stamp

  #-----------------------------------------------------------------------------
  # specify output file name

  output_file = f'{opts.output_root}/atm_srf_{opts.dst_grid}_{cdate}.nc'

  #-----------------------------------------------------------------------------
  # open input files as datasets
  
  ds_map = xr.open_dataset(opts.map_file)
  ds_veg = xr.open_dataset(opts.vegetation_file)
  ds_slw = xr.open_dataset(opts.soil_water_file)

  #-----------------------------------------------------------------------------
  # check that dimension sizes are consistent

  nxy_veg = len(ds_veg.lat) * len(ds_veg.lon)
  nxy_slw = len(ds_slw.lat) * len(ds_slw.lon)
  n_a = len(ds_map.n_a)
  n_b = len(ds_map.n_b)

  if nxy_veg!=nxy_slw or n_a!=nxy_veg or n_a!=nxy_slw:
    err_msg = f'{clr.RED}Input file dimension sizes are inconsistent:{clr.END}'
    err_msg += f' n_a: {n_a}  nxy_veg: {nxy_veg}  nxy_slw: {nxy_slw}'
    raise ValueError(err_msg)

  ntime = 12 # expected length of time dimension

  if len(ds_veg.time)!=ntime:
    err_msg = f'vegetation_file time record error:'
    err_msg += f' expected length = {ntime}, actual length = {len(ds_veg.time)}'
    raise ValueError(err_msg)

  if len(ds_slw.time)!=ntime:
    err_msg = f'soil_water_file time record error:'
    err_msg += f' expected length = {ntime}, actual length = {len(ds_slw.time)}'
    raise ValueError(err_msg)

  if len(ds_veg.lat)!=len(ds_slw.lat):
    err_msg = f'inconsistent latitude coordinate lengths:'
    err_msg += f'  ds_veg.lat = {len(ds_veg.lat)}'
    err_msg += f'  ds_slw.lat = {len(ds_slw.lat)}'
    raise ValueError(err_msg)

  if len(ds_veg.lon)!=len(ds_slw.lon):
    err_msg = f'inconsistent latitude coordinate lengths:'
    err_msg += f'  ds_veg.lon = {len(ds_veg.lon)}'
    err_msg += f'  ds_slw.lon = {len(ds_slw.lon)}'
    raise ValueError(err_msg)

  #-----------------------------------------------------------------------------
  # print some informative stuff

  print(f'''
  File and parameter values:
    {clr.GREEN}map_file        {clr.END}: {opts.map_file}
    {clr.GREEN}vegetation_file {clr.END}: {opts.vegetation_file}
    {clr.GREEN}soil_water_file {clr.END}: {opts.soil_water_file}
    {clr.GREEN}output_root     {clr.END}: {opts.output_root} 
    {clr.GREEN}output_file     {clr.END}: {output_file}
    {clr.GREEN}dst_grid        {clr.END}: {opts.dst_grid}
    {clr.GREEN}map_file n_a    {clr.END}: {n_a}
    {clr.GREEN}map_file n_b    {clr.END}: {n_b}
  ''')

  #-----------------------------------------------------------------------------
  # Load input data
  # (also convert percentage to fraction [0,1])
  LANDMASK_in    = ds_veg['LANDMASK'   ].values
  PCT_LAKE_in    = ds_veg['PCT_LAKE'   ].values/1e2
  PCT_URBAN_in   = ds_veg['PCT_URBAN'  ].values/1e2
  PCT_WETLAND_in = ds_veg['PCT_WETLAND'].values/1e2
  PCT_PFT_in     = ds_veg['PCT_PFT'    ].values/1e2
  SOILW_in       = ds_slw['SOILW'      ].values
  #-----------------------------------------------------------------------------
  # Adjust input data based on land fraction (for consistency with fortran tool)
  nlat = len(ds_veg.lat)
  nlon = len(ds_veg.lon)
  adjust_inputs(ntime, nlat, nlon, LANDMASK_in, PCT_LAKE_in, PCT_URBAN_in, PCT_WETLAND_in, SOILW_in)
  #-----------------------------------------------------------------------------
  # Remap data to the atmosphere grid

  n_s = len(ds_map.n_s)
  wgt = ds_map.S.values
  row = ds_map.row.values-1
  col = ds_map.col.values-1

  shp_1D = (n_b)
  shp_2D = (ntime,n_b)

  SOILW = apply_map_2D( np.zeros(shp_2D), ntime, n_s, wgt, row, col, np.reshape(SOILW_in,(ntime,-1)) )
  
  PCT_LAKE    = apply_map_1D( np.zeros(shp_1D), n_s, wgt, row, col, np.reshape(PCT_LAKE_in   ,-1) )
  PCT_URBAN   = apply_map_1D( np.zeros(shp_1D), n_s, wgt, row, col, np.reshape(PCT_URBAN_in  ,-1) )
  PCT_WETLAND = apply_map_1D( np.zeros(shp_1D), n_s, wgt, row, col, np.reshape(PCT_WETLAND_in,-1) )

  npft = len(ds_veg.PCT_PFT[:,0,0])
  PCT_PFT = np.zeros((npft,n_b))
  for p in range(npft):
    PCT_PFT[p,:] = apply_map_1D( np.zeros(shp_1D), n_s, wgt, row, col, np.reshape(PCT_PFT_in[p,:,:],-1) )

  #-----------------------------------------------------------------------------
  # Compute fields to output

  nclass_landuse = 11
  fraction_landuse = np.zeros([nclass_landuse,n_b])
  total_soilw      = np.zeros([ntime,n_b])

  for i in range(n_b):
    # Calculate total_land as sum as all surface type fractions
    total_land = PCT_LAKE[i] + PCT_URBAN[i] + PCT_WETLAND[i]
    for p in range(npft): total_land += PCT_PFT[p,i]
    # Adjust lake area fraction
    if total_land<1.0:
      PCT_LAKE[i] = PCT_LAKE[i] + (1.0 - total_land)
    # Calculate total_soilw
    fraction_soilw = total_land - ( PCT_LAKE[i] + PCT_WETLAND[i] )
    for t in range(ntime):
      total_soilw[t,i] = total_soilw[t,i] + SOILW[t,i] * fraction_soilw

    # Calculate fraction_landuse - n-1 indexing is used to correspond to fortran indexing
    fraction_landuse[ 1-1,i] = PCT_URBAN[i]
    fraction_landuse[ 2-1,i] = np.sum( PCT_PFT[[n-1 for n in [16,17]]    ,i], axis=0)
    fraction_landuse[ 3-1,i] = np.sum( PCT_PFT[[n-1 for n in [13,14,15]] ,i], axis=0)
    fraction_landuse[ 4-1,i] = np.sum( PCT_PFT[[n-1 for n in [5,6,7,8,9]],i], axis=0)
    fraction_landuse[ 5-1,i] = np.sum( PCT_PFT[[n-1 for n in [2,3,4]]    ,i], axis=0)
    fraction_landuse[ 6-1,i] = PCT_WETLAND[i]
    fraction_landuse[ 7-1,i] = PCT_LAKE[i]
    fraction_landuse[ 8-1,i] = PCT_PFT[1-1,i]
    fraction_landuse[11-1,i] = np.sum( PCT_PFT[[n-1 for n in [10,11,12]] ,i], axis=0)

    # Normalize fraction_landuse if it does not sum to 1 +/- tolerance
    fraction_landuse_tolerance = 0.001
    fraction_landuse_sum = np.sum(fraction_landuse[:,i])
    fraction_landuse_error = fraction_landuse_sum - 1.0
    if np.absolute(fraction_landuse_error) > fraction_landuse_tolerance:
      for c in range(nclass_landuse):
        fraction_landuse[c,i] = fraction_landuse[c,i] / fraction_landuse_sum

  #-----------------------------------------------------------------------------
  # Create output file and add fields

  ds_out = xr.Dataset()
  ds_out['fraction_landuse'] = xr.DataArray(fraction_landuse,dims=['class','ncol'])
  ds_out['soilw']            = xr.DataArray(SOILW           ,dims=['month','ncol'])
  ds_out['lon']              = ds_map.xc_b.rename({'n_b':'ncol'})
  ds_out['lat']              = ds_map.yc_b.rename({'n_b':'ncol'})
  ds_out.to_netcdf(path=output_file,mode='w')

  #-----------------------------------------------------------------------------
  # convert to 64bit_data

  cmd = f'ncks -O --fl_fmt=64bit_data {output_file} {output_file} '
  sp.check_call(cmd,shell=True)

  #-----------------------------------------------------------------------------
  # print status message

  print(f'\nsuccessfully created atmsrf file: {clr.MAGENTA}{output_file}{clr.END}')
  print()

#---------------------------------------------------------------------------------------------------
@numba.njit()
def adjust_inputs(ntime, nlat, nlon, LANDMASK, PCT_LAKE, PCT_URBAN, PCT_WETLAND, SOILW):
  for j in range(nlat):
    for i in range(nlon):
      if LANDMASK[j,i]==0:
        PCT_LAKE[j,i]    = 1.0
        PCT_URBAN[j,i]   = 0.0
        PCT_WETLAND[j,i] = 0.0
        for t in range(ntime): SOILW[t,j,i] = 0.0

@numba.njit()
def apply_map_1D( data_out, n_s, wgt, row, col, data_in ):
  for k in range(n_s):
    data_out[row[k]] = data_out[row[k]] + data_in[col[k]] * wgt[k]
  return data_out

@numba.njit()
def apply_map_2D( data_out, ntime, n_s, wgt, row, col, data_in ):
  for t in range(ntime):
    for k in range(n_s):
      data_out[t,row[k]] = data_out[t,row[k]] + data_in[t,col[k]] * wgt[k]
  return data_out
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  main()
#---------------------------------------------------------------------------------------------------