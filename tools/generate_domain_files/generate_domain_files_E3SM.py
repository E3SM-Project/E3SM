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
output_netcdf_type = 'NETCDF3_64BIT_DATA'
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

  For "tri-grid" configurations of E3SM (land grid is different from atmos/ocn):
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
                  default=1e-3,
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
parser.add_option('--fix-pole',
                  dest='fix_pole',
                  default=False,
                  action='store_true',
                  help='Ensure that pole latitude values are exactly 90 degrees')
(opts, args) = parser.parse_args()
#---------------------------------------------------------------------------------------------------
def main():
  global domain_a_grid_file, domain_b_grid_file, ocn_grid_file, atm_grid_file
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
    cdate = datetime.datetime.now().strftime('%Y%m%d')
  else:
    cdate = opts.date_stamp

  #-------------------------------------------------------------------------------
  # specify output file names

  domain_file_ocn_on_ocn = f'{opts.output_root}/domain.ocn.{opts.ocn_grid}.{cdate}.nc'
  domain_file_lnd_on_atm = f'{opts.output_root}/domain.lnd.{opts.lnd_grid}_{opts.ocn_grid}.{cdate}.nc'
  domain_file_ocn_on_atm = f'{opts.output_root}/domain.ocn.{opts.lnd_grid}_{opts.ocn_grid}.{cdate}.nc'

  # cosmetic clean up of file names
  domain_file_ocn_on_ocn = domain_file_ocn_on_ocn.replace('//','/')
  domain_file_lnd_on_atm = domain_file_lnd_on_atm.replace('//','/')
  domain_file_ocn_on_atm = domain_file_ocn_on_atm.replace('//','/')

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
    {clr.GREEN}fix_pole        {clr.END}: {opts.fix_pole}
  ''')

  #-----------------------------------------------------------------------------
  # open map file as dataset
  
  ds = xr.open_dataset(opts.map_file)

  #-----------------------------------------------------------------------------
  # read grid meta-data from map file

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

  src_grid_rank = len(ds.src_grid_rank.values)
  dst_grid_rank = len(ds.dst_grid_rank.values)
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
    {clr.CYAN}src_grid_rank        {clr.END}: {src_grid_rank}
    {clr.CYAN}dst_grid_rank        {clr.END}: {dst_grid_rank}
    {clr.CYAN}src_grid_dims        {clr.END}: {ds.src_grid_dims.values}
    {clr.CYAN}dst_grid_dims        {clr.END}: {ds.dst_grid_dims.values}
  ''')

  #-----------------------------------------------------------------------------
  # deal with pole points

  if opts.fix_pole and src_grid_rank==1:
    print('  NOTE: the --fix-pole option is not appropriate when the source grid is unstructured (i.e. rank=1)')

  if opts.fix_pole and dst_grid_rank==1:
    print('  NOTE: the --fix-pole option is not appropriate when the destination grid is unstructured (i.e. rank=1)')

  if opts.fix_pole and src_grid_rank==2:
    ni = ds.src_grid_dims[0].values
    nj = ds.src_grid_dims[1].values
    yc_a_tmp = ds['yc_a'].values.reshape([nj,ni])
    yc_a_tmp[ 0,:] = -90.0000000000000000000
    yc_a_tmp[-1,:] =  90.0000000000000000000
    ds['yc_a'][:] = yc_a_tmp.reshape([nj*ni])

  if opts.fix_pole and dst_grid_rank==2:
    ni = ds.dst_grid_dims[0].values
    nj = ds.dst_grid_dims[1].values
    yc_b_tmp = ds['yc_b'].values.reshape([nj,ni])
    yc_b_tmp[ 0,:] = -90.0000000000000000000
    yc_b_tmp[-1,:] =  90.0000000000000000000
    ds['yc_b'][:] = yc_b_tmp.reshape([nj*ni])

  #-----------------------------------------------------------------------------
  # Create ocean domain on ocean grid (domain_file_ocn_on_ocn)

  # Get ocn mask on ocn grid
  omask = get_mask(ds,opts,suffix='_a')
  ofrac = xr.zeros_like(ds['area_a'])

  ds_out = xr.Dataset()

  if src_grid_rank==1:
    ds_out['xc']   = ds['xc_a']  .expand_dims(dim='nj').rename({'n_a':'ni'})
    ds_out['yc']   = ds['yc_a']  .expand_dims(dim='nj').rename({'n_a':'ni'})
    ds_out['xv']   = ds['xv_a']  .expand_dims(dim='nj').rename({'n_a':'ni','nv_a':'nv'})
    ds_out['yv']   = ds['yv_a']  .expand_dims(dim='nj').rename({'n_a':'ni','nv_a':'nv'})
    ds_out['area'] = ds['area_a'].expand_dims(dim='nj').rename({'n_a':'ni'})
    ds_out['frac'] = ofrac       .expand_dims(dim='nj').rename({'n_a':'ni'})
    ds_out['mask'] = omask       .expand_dims(dim='nj').rename({'n_a':'ni'})

  if src_grid_rank==2:
    ni = ds.src_grid_dims[0].values
    nj = ds.src_grid_dims[1].values
    ds_out['xc']   = xr.DataArray( ds['xc_a']  .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['yc']   = xr.DataArray( ds['yc_a']  .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['xv']   = xr.DataArray( ds['xv_a']  .values.reshape([nj,ni,4]), dims=['nj','ni','nv'])
    ds_out['yv']   = xr.DataArray( ds['yv_a']  .values.reshape([nj,ni,4]), dims=['nj','ni','nv'])
    ds_out['area'] = xr.DataArray( ds['area_a'].values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['frac'] = xr.DataArray( ofrac       .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['mask'] = xr.DataArray( omask       .values.reshape([nj,ni]),   dims=['nj','ni'])

  ds_out.to_netcdf(path=domain_file_ocn_on_ocn, mode='w', format=output_netcdf_type)

  print(f'successfully created domain file: {clr.MAGENTA}{domain_file_ocn_on_ocn}{clr.END}')

  #-----------------------------------------------------------------------------
  # Create land and ocean domains on atmosphere grid
  xc   = ds['xc_b']
  yc   = ds['yc_b']
  xv   = ds['xv_b']
  yv   = ds['yv_b']
  area = ds['area_b']

  #-----------------------------------------------------------------------------
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
  omask = xr.ones_like(ds['area_b'],dtype=np.int32)
  lmask = xr.zeros_like(ds['area_b'],dtype=np.int32)
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
    lfrac clipped below / above : {opts.fminval      :+10.6e} / {opts.fmaxval      :+10.6f}
    original lfrac min/max      : {lfrac_min         :+10.6e} / {lfrac_max         :+10.6f}
    final    lfrac min/max      : {lfrac.min().values:+10.6e} / {lfrac.max().values:+10.6f}
    ------------------------------------------------------------------------------------------------
  ''')

  #-----------------------------------------------------------------------------
  # Write land domain on atmosphere grid (domain_file_lnd_on_atm)

  ds_out = xr.Dataset()

  if dst_grid_rank==1:
    ds_out['xc']   = ds['xc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['yc']   = ds['yc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['xv']   = ds['xv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
    ds_out['yv']   = ds['yv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
    ds_out['area'] = ds['area_b'].expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['frac'] = lfrac       .expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['mask'] = lmask       .expand_dims(dim='nj').rename({'n_b':'ni'})

  if dst_grid_rank==2:
    ni = ds.dst_grid_dims[0].values
    nj = ds.dst_grid_dims[1].values
    ds_out['xc']   = xr.DataArray( ds['xc_b']  .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['yc']   = xr.DataArray( ds['yc_b']  .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['xv']   = xr.DataArray( ds['xv_b']  .values.reshape([nj,ni,4]), dims=['nj','ni','nv'])
    ds_out['yv']   = xr.DataArray( ds['yv_b']  .values.reshape([nj,ni,4]), dims=['nj','ni','nv'])
    ds_out['area'] = xr.DataArray( ds['area_b'].values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['frac'] = xr.DataArray( lfrac       .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['mask'] = xr.DataArray( lmask       .values.reshape([nj,ni]),   dims=['nj','ni'])

  add_metadata(ds_out)

  ds_out.to_netcdf(path=domain_file_lnd_on_atm, mode='w', format=output_netcdf_type)

  print(f'successfully created domain file: {clr.MAGENTA}{domain_file_lnd_on_atm}{clr.END}')

  #-----------------------------------------------------------------------------
  # Write ocean domain on atmosphere grid (domain_file_ocn_on_atm)

  ds_out = xr.Dataset()

  if dst_grid_rank==1:
    ds_out['xc']   = ds['xc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['yc']   = ds['yc_b']  .expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['xv']   = ds['xv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
    ds_out['yv']   = ds['yv_b']  .expand_dims(dim='nj').rename({'n_b':'ni','nv_b':'nv'})
    ds_out['area'] = ds['area_b'].expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['frac'] = ofrac       .expand_dims(dim='nj').rename({'n_b':'ni'})
    ds_out['mask'] = omask       .expand_dims(dim='nj').rename({'n_b':'ni'})

  if dst_grid_rank==2:
    ni = ds.dst_grid_dims[0].values
    nj = ds.dst_grid_dims[1].values
    ds_out['xc']   = xr.DataArray( ds['xc_b']  .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['yc']   = xr.DataArray( ds['yc_b']  .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['xv']   = xr.DataArray( ds['xv_b']  .values.reshape([nj,ni,4]), dims=['nj','ni','nv'])
    ds_out['yv']   = xr.DataArray( ds['yv_b']  .values.reshape([nj,ni,4]), dims=['nj','ni','nv'])
    ds_out['area'] = xr.DataArray( ds['area_b'].values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['frac'] = xr.DataArray( ofrac       .values.reshape([nj,ni]),   dims=['nj','ni'])
    ds_out['mask'] = xr.DataArray( omask       .values.reshape([nj,ni]),   dims=['nj','ni'])

  add_metadata(ds_out)

  ds_out.to_netcdf(path=domain_file_ocn_on_atm, mode='w', format=output_netcdf_type)

  print(f'successfully created domain file: {clr.MAGENTA}{domain_file_ocn_on_atm}{clr.END}')
  print()

#---------------------------------------------------------------------------------------------------
def add_metadata(ds):
  global opts, domain_a_grid_file, domain_b_grid_file, ocn_grid_file, atm_grid_file
  # add variable attributes
  ds['xc'] = ds['xc'].assign_attrs({'long_name':'longitude of grid cell center'})
  ds['xc'] = ds['xc'].assign_attrs({'units':'degrees_east'})
  ds['xc'] = ds['xc'].assign_attrs({'bounds':'xv'})
  
  ds['yc'] = ds['yc'].assign_attrs({'long_name':'latitude of grid cell center'})
  ds['yc'] = ds['yc'].assign_attrs({'units':'degrees_north'})
  ds['yc'] = ds['yc'].assign_attrs({'bounds':'yv'})

  ds['xv'] = ds['xv'].assign_attrs({'long_name':'longitude of grid cell verticies'})
  ds['xv'] = ds['xv'].assign_attrs({'units':'degrees_east'})

  ds['yv'] = ds['yv'].assign_attrs({'long_name':'latitude of grid cell verticies'})
  ds['yv'] = ds['yv'].assign_attrs({'units':'degrees_north'})

  ds['area'] = ds['area'].assign_attrs({'long_name':'area of grid cell in radians squared'})
  ds['area'] = ds['area'].assign_attrs({'units':'radian2'})
  ds['area'] = ds['area'].assign_attrs({'coordinates':'xc yc'})

  ds['frac'] = ds['frac'].assign_attrs({'long_name':'fraction of grid cell that is active'})
  ds['frac'] = ds['frac'].assign_attrs({'units':'unitless'})
  ds['frac'] = ds['frac'].assign_attrs({'coordinates':'xc yc'})
  ds['frac'] = ds['frac'].assign_attrs({'filter':f'limit frac to [fminval,fmaxval]; fminval={opts.fminval} fmaxval={opts.fmaxval}'})

  ds['mask'] = ds['mask'].assign_attrs({'long_name':'domain mask'})
  ds['mask'] = ds['mask'].assign_attrs({'units':'unitless'})
  ds['mask'] = ds['mask'].assign_attrs({'coordinates':'xc yc'})
  ds['mask'] = ds['mask'].assign_attrs({'comment':'0 value indicates cell is not active'})

  # add global attributes
  ds.attrs['title']             = 'E3SM domain data'
  ds.attrs['Conventions']       = 'CF-1.0'
  ds.attrs['source_code']       = source_code_meta
  ds.attrs['hostname']          = str(host)
  ds.attrs['history']           = f'created by {user}, '+datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
  ds.attrs['source']            = opts.map_file
  ds.attrs['map_domain_a']      = domain_a_grid_file
  ds.attrs['map_domain_b']      = domain_b_grid_file
  ds.attrs['map_grid_file_ocn'] = ocn_grid_file
  ds.attrs['map_grid_file_atm'] = atm_grid_file
  if opts.fix_pole:
    ds.attrs['filter1']         = '--fix-pole option invoked, yc = -+90 at j=1,j=nj'
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
domain_a_grid_file = 'unknown'
domain_b_grid_file = 'unknown'
ocn_grid_file      = 'unknown'
atm_grid_file      = 'unknown'
#---------------------------------------------------------------------------------------------------
if __name__ == '__main__':
  main()
#---------------------------------------------------------------------------------------------------
