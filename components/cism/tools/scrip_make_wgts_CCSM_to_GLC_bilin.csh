#!/bin/csh -fv 
#===============================================================================
# This scrip will create standard CCSM area-average and bilinear mapping data 
# files in the directory in which this script is run.
# 
# This script assumes
# 1) that the scrip source code has already been compiled appropriately
# 2) that both source & destination grid definition files already exist
#===============================================================================

set echo off
set verbose

#-------------------------------------------------------------------------------
# edit these csh vars to define what mapping files will be created
#-------------------------------------------------------------------------------

set SCRIPDIR   = /blhome/jwolfe/scrip/scrip1.4_pop-roms_01
set GRIDDIR    = /fis/cgd/cseg/csm/mapping/grids

set gridFn_glc = gland10km_scrip.nc
set gridFn_lnd = $GRIDDIR/T31_040122.nc

set DATE = "`date +%y%m%d`"

set mapFn_l2g_bilin  = map_T31_to_GLC_Gland10km_bilin_da_$DATE.nc
set map1_name        = 'T31 to GLC test (unmasked) mapping'
set map_method_bilin = bilinear
set normalize_opt    = destarea

#===============================================================================
# do not edit below this point
#===============================================================================

#-------------------------------------------------------------------------------
# create scrip input namelist files
#-------------------------------------------------------------------------------
# notes: 
# o scrip will only recognize a domain mask for "grid1", hence grid1 = ocn grid
# o one invocation of scrip creates two maps...
#   map1: grid1 -> grid2 
#   map2: grid2 -> grid1
#-------------------------------------------------------------------------------

cat >! scrip_in_bilin << EOF
&remap_inputs
    num_maps        = 1
    grid1_file      = '$gridFn_lnd'
    grid2_file      = '$gridFn_glc'
    interp_file1    = '$mapFn_l2g_bilin'
    map1_name       = '$map1_name'
    map_method      = '$map_method_bilin'
    normalize_opt   = '$normalize_opt'
    output_opt      = 'ncar-csm'
    restrict_type   = 'latlon'
    num_srch_bins   = 10
    luse_grid1_area = .false.
    luse_grid2_area = .false.
    grid2_periodic  = .false.
/
EOF

#-------------------------------------------------------------------------------
# run scrip: create the netCDF mapping data files
#-------------------------------------------------------------------------------
# notes: 
# o scrip namelist file must be named "scrip_in"
#-------------------------------------------------------------------------------

cp scrip_in_bilin scrip_in
$SCRIPDIR/scrip  

#-------------------------------------------------------------------------------
# append documentation to the mapping data files (netCDF global attribute)
#-------------------------------------------------------------------------------

foreach FILE ( $mapFn_l2g_bilin)
  ncatted -O -h -a Created_by,global,a,c,"`whoami`, `date`" $FILE
  ncatted -O -h -a 1D_grid_indexing,global,a,c,"if n is 1D index, i runs fast, j runs slow: n=(j-1)*fast_grid_dim+i" $FILE
  ncatted -O -h -a grid_file_glc,global,a,c,"$gridFn_glc" $FILE
  ncatted -O -h -a grid_file_ldn,global,a,c,"$gridFn_lnd" $FILE
end


exit(0)
