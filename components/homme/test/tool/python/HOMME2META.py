#!/usr/bin/env python3
import os
import xarray as xr
import datetime
user, host = os.getenv('USER'), os.getenv('HOST')
source_code_meta = 'HOMME2META.py'
output_netcdf_type = 'NETCDF3_64BIT_DATA'
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
verbose_indent = ' '*2
#---------------------------------------------------------------------------------------------------
usage = '''
python HOMME2META.py --src_file <src_file> --dst_file <dst_file> 

Purpose:
  This script reads a HOMME grid template file and writes out a "latlon" format grid 
  description file of the np4/GLL grid.  The output file contains a list of all the 
  unique GLL nodes (no duplicate degrees of freedom) and their subcell connectivity 
  (dividing each spectral element into a (np-1) x (np-1) subcells with GLL nodes as their corners)
  
  This "latlon" template file can be used by plotting programs to plot GLL vertex
  data from the dycore (as opposed to the PG2 finite volume cell centered data from the physics)
  It is also used by the MOAB disk averaging remapping algorithm.

Environment
  
  This requires libraries such as xarray, which is included in the E3SM unified environment:
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
                  help='Output latlon grid file')
(opts, args) = parser.parse_args()


def main():
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
  ds_in = xr.open_dataset(opts.src_file)

  
  # Remove existing output if it exists
  if os.path.exists(opts.dst_file):
        os.remove(opts.dst_file)


  # Compute grid size
  lat = ds_in["lat"]
  print(f"grid_size = {lat.size}")

  # Print lon/lat min/max
  lon = ds_in["lon"]
  print(f"lon min/max = {float(lon.min())} {float(lon.max())}")
  print(f"lat min/max = {float(lat.min())} {float(lat.max())}")

  # Check for corners variable and slice first 4 rows
  print("reading corners...")
  corners_var = ds_in["corners"].values[0:4, :]   # only first 4 rows
  corners_int = corners_var.astype(int)
  
  # Prepare output dataset
  print("opening dataset...")
  ds_out = xr.Dataset()

  # Copy variables lat, lon, area
  ds_out["area"]  = ds_in["area"]
  ds_out["lat"]  = ds_in["lat"]  # (("ncol"),ds_in["lat"])
  ds_out["lon"]  = ds_in["lon"]  # (("ncol"),ds_in["lon"].values)
  # element_corners has 4 corners per cell; use "ncorners" and "ncells" as dimension names
  ds_out["element_corners"] = (("ncorners","ncells"), corners_int)

  # Global attributes
  ds_out.attrs = ds_in.attrs.copy()
  ds_out.attrs['title']             = 'HOMME generated GLL grid data, latlon format'
  ds_out.attrs['hostname']          = str(host)
  ds_out.attrs['history']           = f'created by {user}, '+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  ds_out.attrs['src_file']          = opts.src_file
  ds_out.attrs['source_code']       = source_code_meta

  # Write to netCDF using the specified format
  print("writing dataset...")  
  ds_out.to_netcdf(opts.dst_file, format=output_netcdf_type)
  print(f"Wrote output file: {opts.dst_file} with format {output_netcdf_type}")

if __name__ == "__main__":
    main()
