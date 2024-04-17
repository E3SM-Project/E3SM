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
The map file used to generate the domain files can be created a few different ways.
For a typical E3SM configuration we recommend using a conservative, monotone map.
Here is an example command that can be used to generate one as of NCO version 5.2.2

  ncremap -5 -a traave --src_grd=${OCN_GRID} --dst_grd=${ATM_GRID} --map_file=${MAP_FILE}

'''
#---------------------------------------------------------------------------------------------------
import datetime, os, numpy as np, xarray as xr, numba, itertools
user, host = os.getenv('USER'), os.getenv('HOST')
source_code_meta = 'generate_domain_E3SM.py'
#---------------------------------------------------------------------------------------------------
class clr:END,RED,GREEN,MAGENTA,CYAN = '\033[0m','\033[31m','\033[32m','\033[35m','\033[36m'
#---------------------------------------------------------------------------------------------------
usage = '''
python generate_domain_files_E3SM.py  -m <map_file> 
                                      -o <ocn_grid_name> 
                                      -l <lnd_grid_name>
                                      [--output-root <path>]
                                      [--date-stamp <date string>]
                                      [--fminval <fminval>]
                                      [--fmaxval <fmaxval>]
                                      [--set-omask]

Purpose:
  For "bi-grid" configurations of E3SM (land grid is same as atmos):
  Given a mapping file from the ocean grid (where the mask is defined) 
  to the atmosphere grid, this tool creates land and ocean domain files 
  needed by data model components (ex. datm, dlnd, docn)

  For "tr-grid" configurations of E3SM (land grid is different from atmos/ocn):
  In addition to running this tool with the ocn->atm map as above,
  a second iteration is needed with a similar ocn->lnd map.

Environment
  
  This tool requires a few special packages, such as xarray, numba, and itertools.
  These are all included in the E3SM unified environment:
  https://e3sm.org/resources/tools/other-tools/e3sm-unified-environment/

  Otherwise a simple conda environment can be created:
  conda create --name example_env --channel conda-forge xarray numpy numba scikit-learn netcdf4

The following output domain files are created:

  domain.lnd.<gridlnd>_<gridocn>.<date_stamp>.nc
    land domain file on the land/atmos grid with a land fraction 
    corresponding to (1-ocnfrac) mask mapped to the land grid

  domain.ocn.<gridlnd>_<gridocn>.<date_stamp>.nc
    ocean domain on the land/atmos grid with an ocean fraction based 
    on the ocean grid mask mapped to the land/atmos grid for when 
    atm,lnd,ice,ocn are all on the same grid 
    (not compatible with MPAS sea-ice)

  domain.ocn.<gridocn>.<date_stamp>.nc
    ocean domain on the ocean grid 
'''
from optparse import OptionParser
parser = OptionParser(usage=usage)
parser.add_option('-m',
                  dest='map_file',
                  default=None,
                  help='Input mapping file from ocean to atmosphere grid - '+\
                       'This is the primary source from which the ocean and'+\
                       'land domains will be determined')
parser.add_option('-l',
                  dest='lnd_grid',
                  default=None,
                  help='Output land grid name')
parser.add_option('-o',
                  dest='ocn_grid',
                  default=None,
                  help='Output ocean grid name')
parser.add_option('--output-root',
                  dest='output_root',
                  default='./',
                  help='Output path for domain files')
parser.add_option('--date-stamp',
                  dest='date_stamp',
                  default=None,
                  help='Creation date stamp for domain files')
parser.add_option('--fminval',
                  dest='fminval',
                  default=0, 
                  help='Minimum allowable land fraction (reset to 0 below fminval)')
parser.add_option('--fmaxval',
                  dest='fmaxval',
                  default=1,
                  help='Maximum allowable land fraction (reset to 1 above fmaxval)')
parser.add_option('--set-omask',
                  dest='set_omask',
                  default=False,
                  action='store_true',
                  help='If True then an ocean mask is not required and will '+\
                       'simply be set to a vector of 1\'s  if mask_a is not '+\
                       'present in the input mapping file. If --set-omask is '+\
                       'omitted, then mask_a is required and an error will be '+\
                       'raised if it does not exist in the input mapping file.')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
def main():

  #-------------------------------------------------------------------------------
  # check for valid input arguments

  if opts.map_file is None: 
    raise ValueError(f'{clr.RED}input map file was not specified{clr.END}')
  if opts.lnd_grid is None: 
    raise ValueError(f'{clr.RED}land grid name was not specified{clr.END}')
  if opts.ocn_grid is None: 
    raise ValueError(f'{clr.RED}ocean grid name was not specified{clr.END}')
  if not os.path.exists(opts.output_root) :
    raise ValueError(f'{clr.RED}Output root path does not exist{clr.END}')

  #-------------------------------------------------------------------------------
  # Set date stamp for file name

  if opts.date_stamp is None:
    cdate = datetime.datetime.utcnow().strftime('%Y%m%d')
  else:
    cdate = opts.date_stamp

  #-------------------------------------------------------------------------------
  # specify output file names

  domain_file_ocn_on_ocn = f'{opts.output_root}/test.domain.ocn.{opts.ocn_grid}.{cdate}.nc'
  domain_file_lnd_on_atm = f'{opts.output_root}/test.domain.lnd.{opts.lnd_grid}_{opts.ocn_grid}.{cdate}.nc'
  domain_file_ocn_on_atm = f'{opts.output_root}/test.domain.ocn.{opts.lnd_grid}_{opts.ocn_grid}.{cdate}.nc'

  #-----------------------------------------------------------------------------
  # print some informative stuff

  print(f'''
  Input files and parameter values:
    {clr.GREEN}map_file        {clr.END}: {opts.map_file}
    {clr.GREEN}lnd_grid        {clr.END}: {opts.lnd_grid}
    {clr.GREEN}ocn_grid        {clr.END}: {opts.ocn_grid}
    {clr.GREEN}fminval         {clr.END}: {opts.fminval}
    {clr.GREEN}fmaxval         {clr.END}: {opts.fmaxval}
    {clr.GREEN}set_omask       {clr.END}: {opts.set_omask}
  ''')

  #-----------------------------------------------------------------------------
  # open map file as dataset
  
  ds = xr.open_dataset(opts.map_file)

  #-----------------------------------------------------------------------------
  # read grid meta-data from map file

  domain_a_grid_file = 'unknown'
  domain_b_grid_file = 'unknown'
  ocn_grid_file      = 'unknown'
  atm_grid_file      = 'unknown'

  if 'domain_a' in ds.attrs.keys(): domain_a_grid_file = ds.attrs['domain_a']
  if 'domain_b' in ds.attrs.keys(): domain_b_grid_file = ds.attrs['domain_b']

  if 'grid_file_ocn' in ds.attrs.keys():
    ocn_grid_file = ds.attrs['grid_file_ocn']
  else:
    ocn_grid_file = ds.attrs['grid_file_src']

  if 'grid_file_ocn' in ds.attrs.keys():
    atm_grid_file = ds.attrs['grid_file_atm']
  else:
    atm_grid_file = ds.attrs['grid_file_dst']

  #-----------------------------------------------------------------------------
  # print some useful information from the map file

  print(f'''
  Grid information from map file:
    {clr.CYAN}domain_a file        {clr.END}: {domain_a_grid_file}
    {clr.CYAN}domain_b file        {clr.END}: {domain_b_grid_file}
    {clr.CYAN}ocn_grid_file        {clr.END}: {ocn_grid_file}
    {clr.CYAN}atm_grid_file        {clr.END}: {atm_grid_file}
    {clr.CYAN}ocn grid size   (n_a){clr.END}: {len(ds.n_a)}
    {clr.CYAN}atm grid size   (n_b){clr.END}: {len(ds.n_b)}
    {clr.CYAN}sparse mat size (n_s){clr.END}: {len(ds.n_s)}
  ''')

  #-----------------------------------------------------------------------------
  # Create ocean domain on ocean grid (domain_file_ocn_on_ocn)

  # Get ocn mask on ocn grid
  omask = get_mask(ds,opts,suffix='_a')
  ofrac = xr.zeros_like(ds['area_a'])

  ds_out = xr.Dataset()
  ds_out['xc']   = ds['xc_a']  .expand_dims(dim='nj').rename({'n_a':'ni'})
  ds_out['yc']   = ds['yc_a']  .expand_dims(dim='nj').rename({'n_a':'ni'})
  ds_out['xv']   = ds['xv_a']  .expand_dims(dim='nj').rename({'n_a':'ni','nv_a':'nv'})
  ds_out['yv']   = ds['yv_a']  .expand_dims(dim='nj').rename({'n_a':'ni','nv_a':'nv'})
  ds_out['area'] = ds['area_a'].expand_dims(dim='nj').rename({'n_a':'ni'})
  ds_out['frac'] = ofrac       .expand_dims(dim='nj').rename({'n_a':'ni'})
  ds_out['mask'] = omask       .expand_dims(dim='nj').rename({'n_a':'ni'})

  add_metadata(ds_out)

  ds_out.to_netcdf(path=domain_file_ocn_on_ocn,mode='w')

  print(f'successfully created domain file: {clr.MAGENTA}{domain_file_ocn_on_ocn}{clr.END}')

  #-----------------------------------------------------------------------------
  # Create land and ocean domains on atmosphere grid

  xc   = ds['xc_b']
  yc   = ds['yc_b']
  xv   = ds['xv_b']
  yv   = ds['yv_b']
  area = ds['area_b']

  mask_a = get_mask(ds,opts,suffix='_a')
  frac_a = xr.where( mask_a!=0, xr.ones_like(ds['area_a']), xr.zeros_like(ds['area_a']) )

  # compute ocn fraction on atm grid
  ofrac = compute_ofrac_on_atm( len(ds['n_s']), np.zeros(ds['area_b'].shape), 
                                frac_a.values, ds['S'].values, 
                                ds['row'].values-1, ds['col'].values-1 )
  ofrac = xr.DataArray(ofrac,dims=['n_b'])

  # lfrac is area frac of mask "_a" on grid "_b" or float(mask)
  lfrac = xr.zeros_like(ds['area_b'])

  # convert to land fraction
  lfrac_min = opts.fmaxval
  lfrac_max = opts.fminval
  omask = xr.ones_like(ds['area_b'],dtype=int32)
  lmask = xr.zeros_like(ds['area_b'],dtype=int32)
  lfrac = 1 - ofrac
  lfrac_min = lfrac.min().values
  lfrac_max = lfrac.max().values
  lfrac = xr.where( lfrac>opts.fmaxval, 1, lfrac )
  lfrac = xr.where( lfrac<opts.fminval, 0, lfrac )
  ofrac = 1 - lfrac
  lmask = xr.where( lfrac!=0, 1, lmask )
  omask = xr.where( ofrac==0, 0, omask )

  print(f'''
    ------------------------------------------------------------------------------------------------
    {clr.RED}IMPORTANT{clr.END}: note original min/max frac and decide if acceptable
    lfrac clipped below / above : {opts.fminval      :+8.4e} / {opts.fmaxval      :+8.4f}
    original lfrac min/max      : {lfrac_min         :+8.4e} / {lfrac_max         :+8.4f}
    final    lfrac min/max      : {lfrac.min().values:+8.4e} / {lfrac.max().values:+8.4f}
    ------------------------------------------------------------------------------------------------
  ''')

  #-----------------------------------------------------------------------------
  # Write land domain on atmosphere grid (domain_file_lnd_on_atm)

  ds_out = xr.Dataset()
  ds_out['xc']   = ds['xc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['yc']   = ds['yc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['xv']   = ds['xv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
  ds_out['yv']   = ds['yv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
  ds_out['area'] = ds['area_b'].expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['frac'] = lfrac       .expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['mask'] = lmask       .expand_dims(dim='nj').rename({'n_b':'ni'})

  add_metadata(ds_out)

  ds_out.to_netcdf(path=domain_file_lnd_on_atm,mode='w')

  print(f'successfully created domain file: {clr.MAGENTA}{domain_file_lnd_on_atm}{clr.END}')

  #-----------------------------------------------------------------------------
  # Write ocean domain on atmosphere grid (domain_file_ocn_on_atm)

  ds_out = xr.Dataset()
  ds_out['xc']   = ds['xc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['yc']   = ds['yc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['xv']   = ds['xv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
  ds_out['yv']   = ds['yv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
  ds_out['area'] = ds['area_b'].expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['frac'] = ofrac       .expand_dims(dim='nj').rename({'n_b':'ni'})
  ds_out['mask'] = omask       .expand_dims(dim='nj').rename({'n_b':'ni'})

  add_metadata(ds_out)

  ds_out.to_netcdf(path=domain_file_ocn_on_atm,mode='w')

  print(f'successfully created domain file: {clr.MAGENTA}{domain_file_ocn_on_atm}{clr.END}')
  print()

#---------------------------------------------------------------------------------------------------
def add_metadata(ds):
  # add variable attirbutes
  ds_out['xc'] = ds_out['xc'].assign_attrs({'long_name':'longitude of grid cell center'})
  ds_out['xc'] = ds_out['xc'].assign_attrs({'units':'degrees_east'})
  ds_out['xc'] = ds_out['xc'].assign_attrs({'bounds':'xv'})
  
  ds_out['yc'] = ds_out['yc'].assign_attrs({'long_name':'latitude of grid cell center'})
  ds_out['yc'] = ds_out['yc'].assign_attrs({'units':'degrees_north'})
  ds_out['yc'] = ds_out['yc'].assign_attrs({'bounds':'yv'})

  ds_out['xv'] = ds_out['xv'].assign_attrs({'long_name':'longitude of grid cell verticies'})
  ds_out['xv'] = ds_out['xv'].assign_attrs({'units':'degrees_east'})

  ds_out['yv'] = ds_out['yv'].assign_attrs({'long_name':'latitude of grid cell verticies'})
  ds_out['yv'] = ds_out['yv'].assign_attrs({'units':'degrees_north'})

  ds_out['area'] = ds_out['area'].assign_attrs({'long_name':'area of grid cell in radians squared'})
  ds_out['area'] = ds_out['area'].assign_attrs({'units':'radian2'})
  ds_out['area'] = ds_out['area'].assign_attrs({'coordinates':'xc yc'})

  ds_out['frac'] = ds_out['frac'].assign_attrs({'long_name':'fraction of grid cell that is active'})
  ds_out['frac'] = ds_out['frac'].assign_attrs({'units':'unitless'})
  ds_out['frac'] = ds_out['frac'].assign_attrs({'coordinates':'xc yc'})
  ds_out['frac'] = ds_out['frac'].assign_attrs({'filter':f'limit frac to [fminval,fmaxval]; fminval={opts.fminval} fmaxval={opts.fmaxval}'})

  ds_out['mask'] = ds_out['mask'].assign_attrs({'long_name':'domain mask'})
  ds_out['mask'] = ds_out['mask'].assign_attrs({'units':'unitless'})
  ds_out['mask'] = ds_out['mask'].assign_attrs({'coordinates':'xc yc'})
  ds_out['mask'] = ds_out['mask'].assign_attrs({'comment':'0 value indicates cell is not active'})

  # add global attributes
  ds.attrs['title']             = 'E3SM domain data'
  ds.attrs['Conventions']       = 'CF-1.0'
  ds.attrs['source_code']       = source_code_meta
  ds.attrs['hostname']          = host
  ds.attrs['history']           = f'created by {user}, '+datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
  ds.attrs['source']            = opts.map_file
  ds.attrs['map_domain_a']      = domain_a_grid_file
  ds.attrs['map_domain_b']      = domain_b_grid_file
  ds.attrs['map_grid_file_ocn'] = ocn_grid_file
  ds.attrs['map_grid_file_atm'] = atm_grid_file
#---------------------------------------------------------------------------------------------------
def get_mask(ds,opts,suffix):
  mask_var_name = 'mask'+suffix
  if mask_var_name in ds.variables:
    mask_tmp = ds[mask_var_name]
  else:
    if opts.set_omask:
      # set mask=1 everywhere
      mask_tmp = xr.ones_like(ds['area'+suffix])
    else:
      raise ValueError(f'{clr.RED}{mask_var_name} does not exist in input map file{clr.END}')
  return mask_tmp
#---------------------------------------------------------------------------------------------------
@numba.njit()
def compute_ofrac_on_atm( n_s, ofrac, frac_a, S, row, col ):
  for k in range(n_s):
    ofrac[row[k]] = ofrac[row[k]] + frac_a[col[k]] * S[k]
  return ofrac
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  main()
#---------------------------------------------------------------------------------------------------
