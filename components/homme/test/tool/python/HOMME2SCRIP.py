#!/usr/bin/env python3
#---------------------------------------------------------------------------------------------------
'''
This is a replacement for the legacy HOMME2SCRIP.ncl tool created for CESM. 
Most legacy functionality is reproduced

Created May, 2025 by Walter Hannah (LLNL) 
'''
#---------------------------------------------------------------------------------------------------
import datetime, os, numpy as np, xarray as xr, numba
user, host = os.getenv('USER'), os.getenv('HOST')
source_code_meta = 'HOMME2SCRIP.py'
output_netcdf_type = 'NETCDF3_64BIT_DATA'
from time import perf_counter
enable_timers = True
# ------------------------------------------------------------------------------
class tclr: END,GREEN,YELLOW,MAGENTA,CYAN='\033[0m','\033[32m','\033[33m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
verbose_indent = ' '*2
# ------------------------------------------------------------------------------
def print_timer(timer_start,description='',indent=None):
  global verbose_indent
  if indent is None: indent = verbose_indent
  etime = perf_counter()-timer_start
  time_str = f'{etime:10.1f} sec'
  if etime>60: time_str += f' ({(etime/60):4.1f} min)'
  print(f'\n{indent}{tclr.YELLOW}{description:40} elapsed time: {time_str}{tclr.END}')
#---------------------------------------------------------------------------------------------------
usage = '''
python HOMME2SCRIP.py -i <src_file> 
                      -o <dst_file> 
                      --ne <ne>
                      --np <np>

Purpose:
  This script reads a HOMME grid template file and writes out a SCRIP format grid description file of the np4/GLL grid.
  
  HOMME np4 grid template files are produced by a two step procedure, which first requires running homme_tool, and then this script to convert the output into SCRIP format. This procedure is only needed for np4 files due to their use of vertex data. For cell centered pg2 files, one should instead use TempestRemap to create a grid description file. This is particularly useful when remapping topography data with cube_to_target, which can be much faster than remapping with tools like NCO due to the large size of the input topography data.

Environment
  
  This requires libraries such as xarray, which included in the E3SM unified environment:
  https://e3sm.org/resources/tools/other-tools/e3sm-unified-environment/

  Otherwise a simple conda environment can be created:
  conda create --name example_env --channel conda-forge xarray numpy netcdf4

'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('--src_file',
                  dest='src_file',
                  default=None,
                  help='Input HOMME grid template file')
parser.add_option('--dst_file',
                  dest='dst_file',
                  default=None,
                  help='Output scrip grid file')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
def main():
  if enable_timers: timer_start_all = perf_counter()
  #-------------------------------------------------------------------------------
  # check for valid input arguments
  if opts.src_file is None: raise ValueError(f'{clr.RED}src_file argument was not specified{clr.END}')
  if opts.dst_file is None: raise ValueError(f'{clr.RED}dst_file argument was not specified{clr.END}')

  #-----------------------------------------------------------------------------
  # print some informative stuff
  print()
  print(verbose_indent+f'{clr.GREEN}Input arguments:{clr.END}')
  print(verbose_indent+f'  {clr.CYAN}src_file{clr.END}: {opts.src_file}')
  print(verbose_indent+f'  {clr.CYAN}dst_file{clr.END}: {opts.dst_file}')

  #-----------------------------------------------------------------------------
  # open input file as dataset
  ds = xr.open_dataset(opts.src_file,decode_timedelta=True)

  #-----------------------------------------------------------------------------
  # check for variables we need
  if 'lat'    not in ds: raise ValueError(f'{clr.RED}required variable missing from input file:{clr.END} lat')
  if 'lon'    not in ds: raise ValueError(f'{clr.RED}required variable missing from input file:{clr.END} lon')
  if 'area'   not in ds: raise ValueError(f'{clr.RED}required variable missing from input file:{clr.END} area')
  if 'cv_lat' not in ds: raise ValueError(f'{clr.RED}required variable missing from input file:{clr.END} cv_lat')
  if 'cv_lon' not in ds: raise ValueError(f'{clr.RED}required variable missing from input file:{clr.END} cv_lon')

  #-----------------------------------------------------------------------------
  # remove variables we don't need - keep lev since it holds the corner values
  if 'time'    in ds: ds = ds.isel(time=0,drop=True)
  if 'ilev'    in ds: ds = ds.isel(ilev=0,drop=True)
  if 'hyam'    in ds: ds = ds.drop_vars('hyam')
  if 'hybm'    in ds: ds = ds.drop_vars('hybm')
  if 'hyai'    in ds: ds = ds.drop_vars('hyai')
  if 'hybi'    in ds: ds = ds.drop_vars('hybi')
  if 'corners' in ds: ds = ds.drop_vars('corners')

  #-----------------------------------------------------------------------------
  # use the number of valid corner locations to determine max corners across grid
  num_corners = 0
  kmax = ds.cv_lon['lev'].shape[0] # lev value hold corner values

  for k in range(kmax):
    # use max of logitude values to check if any corner values exist for k
    max_lon = np.max( np.absolute( ds.cv_lon.isel(lev=k).values ) )
    # if max_lon is zero then no more corners exist and we can use k as max # of corners
    if ( max_lon<0.000000001 and num_corners==0): num_corners = k

  print()
  print(verbose_indent+f'{clr.GREEN}Unstructured control volumes max number of corners:{clr.END} {num_corners}')

  #-----------------------------------------------------------------------------
  # print min/max of coordinates
  print()
  print(verbose_indent+f'{clr.GREEN}Sanity check for coordinate bounds:{clr.END}')
  print(verbose_indent+f'  lon    min/max: {np.min(ds.lon.values)   :8.4f} / {np.max(ds.lon.values)   :8.4f}')
  print(verbose_indent+f'  lat    min/max: {np.min(ds.lat.values)   :8.4f} / {np.max(ds.lat.values)   :8.4f}')
  print(verbose_indent+f'  cv_lon min/max: {np.min(ds.cv_lon.values):8.4f} / {np.max(ds.cv_lon.values):8.4f}')
  print(verbose_indent+f'  cv_lat min/max: {np.min(ds.cv_lat.values):8.4f} / {np.max(ds.cv_lat.values):8.4f}')

  #-----------------------------------------------------------------------------
  # Create output dataset
  if enable_timers: timer_start = perf_counter()

  ds_out = ds.rename({'ncol'  :'grid_size',\
                      'area'  :'grid_area',\
                      'lev'   :'grid_corners',\
                      'lat'   :'grid_center_lat',\
                      'lon'   :'grid_center_lon',\
                      'cv_lat':'grid_corner_lat',\
                      'cv_lon':'grid_corner_lon',\
                    }).isel(grid_corners=slice(0,num_corners))


  ds_out['grid_area'] = ds_out['grid_area'].assign_attrs(units='radians^2')
  ds_out['grid_area'] = ds_out['grid_area'].assign_attrs(long_name='area weights')

  ds_out['grid_center_lat'] = ds_out['grid_center_lat'].assign_attrs(units='degrees')
  ds_out['grid_center_lon'] = ds_out['grid_center_lon'].assign_attrs(units='degrees')
  ds_out['grid_corner_lat'] = ds_out['grid_corner_lat'].assign_attrs(units='degrees')
  ds_out['grid_corner_lon'] = ds_out['grid_corner_lon'].assign_attrs(units='degrees')

  #-----------------------------------------------------------------------------
  # transpose grid data
  for v in ds_out.variables:
    if 'grid_corners' in ds_out[v].dims:
      ds_out[v] = ds_out[v].transpose('grid_size','grid_corners',missing_dims='ignore')

  # remove the grid_corners variable:
  del ds_out['grid_corners']

  ds_out.load()

  if enable_timers: print_timer(timer_start,description='create output dataset')

  #-----------------------------------------------------------------------------
  @numba.njit()
  def pole_corner_check(ncol,center_lat,center_lon,corner_lat,corner_lon):
    num_modified = 0
    for i in range(ncol):
      pole_dist = np.absolute( 90 - np.absolute( center_lat[i] ) )
      if ( pole_dist < 1e-9 ):
        # 1 2 3 4 -> 1 4 3 2    swap pos 1,3
        tmp_corner_lon = corner_lon[i,:]
        tmp_corner_lat = corner_lat[i,:]
        corner_lon[i,1] = tmp_corner_lon[3]
        corner_lon[i,3] = tmp_corner_lon[1]
        corner_lat[i,1] = tmp_corner_lat[3]
        corner_lat[i,3] = tmp_corner_lat[1]
        num_modified += 1
    return num_modified
  #-----------------------------------------------------------------------------
  # Fix orientation at pole points
  print()
  print(verbose_indent+f'{clr.GREEN}Checking pole coordinates...{clr.END}')

  if enable_timers: timer_start = perf_counter()

  tmp_ncol       = len(ds_out['grid_size'].values)
  tmp_center_lat = ds_out['grid_center_lat'].values
  tmp_center_lon = ds_out['grid_center_lon'].values
  tmp_corner_lat = ds_out['grid_corner_lat'].values
  tmp_corner_lon = ds_out['grid_corner_lon'].values

  num_modified = pole_corner_check(tmp_ncol,tmp_center_lat,tmp_center_lon,tmp_corner_lat,tmp_corner_lon)

  if num_modified>0:
    print(verbose_indent+(' '*2)+f'{clr.GREEN}{num_modified} pole point modified{clr.END}')
  else:
    print(verbose_indent+(' '*2)+f'No pole points modified')

  if enable_timers: print_timer(timer_start,description='pole check')

  #-----------------------------------------------------------------------------
  # add imask and grid_dims to output datasest
  ds_out['grid_imask'] = xr.ones_like(ds_out['grid_size'],dtype=int)
  # HOMME grids are always 1D unstructured:
  ds_out['grid_dims'] = xr.DataArray([1],dims=['grid_rank'])

  #-----------------------------------------------------------------------------
  # add global attributes
  ds_out.attrs['title']             = 'HOMME generated np4 SCRIP grid data'
  ds_out.attrs['Conventions']       = 'CF-1.0'
  ds_out.attrs['source_code']       = source_code_meta
  ds_out.attrs['hostname']          = str(host)
  ds_out.attrs['history']           = f'created by {user}, '+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  ds_out.attrs['src_file']          = opts.src_file

  #-----------------------------------------------------------------------------
  # write grid data out to netcdf file
  print()
  print(verbose_indent+f'{clr.GREEN}Writing output grid data...{clr.END}')

  if enable_timers: timer_start = perf_counter()

  ds_out.to_netcdf(path=opts.dst_file, mode='w', format=output_netcdf_type)

  if enable_timers: print_timer(timer_start,description='write to file')

  #-----------------------------------------------------------------------------
  # final print statements
  print()
  print(verbose_indent+f'{clr.GREEN}Successfully created file:{clr.END} {opts.dst_file}')

  if enable_timers: print_timer(timer_start_all,description='total time')

#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  main()
#---------------------------------------------------------------------------------------------------
