#!/bin/csh -f

#set verbose

set OCN_GRIDS = ( gx3v5 gx3v6 gx1v5 tx1v1 tx0.1v2)
set COMPSETS = (C B)

#--------------------------------------------------------------------------
#  Define environment variables and testing directories
#--------------------------------------------------------------------------

setenv OUTROOT .

setenv rundir  $OUTROOT/testrun_dir
if !(-e $rundir) mkdir $rundir
setenv nmldir  $OUTROOT/nml_dir
if !(-e $nmldir) mkdir $nmldir

setenv INPUT  $OUTROOT/$rundir/input
if !(-e $INPUT) mkdir $INPUT

#--------------------------------------------------------------------------
#  Test by cycling over allowable grids and compsets
#--------------------------------------------------------------------------
foreach grid ($OCN_GRIDS)
foreach compset ($COMPSETS)


if ($compset == 'B') then
 setenv OCN_COUPLING full
 setenv OCN_ICE_FORCING active
else if ($compset == 'C') then
 setenv OCN_COUPLING partial
 setenv OCN_ICE_FORCING inactive
else
 exit -99
endif

setenv OCN_GRID $grid

echo $compset $OCN_GRID

setenv CASE pop2_in_testcase_${OCN_GRID}_${compset}
setenv LID 123456
setenv POP2_NMLFILE $nmldir/pop2_in_${OCN_GRID}_${compset}
setenv NPROCS_CLINIC 4
setenv NPROCS_TROPIC 4
if ($OCN_GRID == tx0.1v2) then
setenv OCN_CDF64 TRUE
setenv OCN_NCPL  4
else
setenv OCN_CDF64 FALSE
setenv OCN_NCPL  1
endif
setenv IYEAR0   0
setenv IMONTH0  1
setenv IDAY0    1
setenv IHOUR0   0
setenv OCN_CHL_TYPE diagnostic

setenv runtype  'ccsm_startup'

setenv output_L  $rundir
setenv output_r  $rundir/$CASE.pop.r
setenv output_h  $rundir/$CASE.pop.h
setenv output_d  $rundir/$CASE.pop.d
setenv log_filename   ocn.log.$LID
setenv pop2_pointer  $rundir/rpointer.ocn
setenv depth_accel_filename        ${grid}_depth_accel
setenv bathymetry_filename         ${grid}_bathymetry
setenv bottom_cell_filename        ${grid}_bottom_cell
setenv chl_filename                ${grid}_chl
setenv horiz_grid_filename         ${grid}_horiz_grid
setenv init_ts_filename            ${grid}_init_ts
setenv region_ids_filename         ${grid}_region_ids
setenv regionmask_filename         ${grid}_regionmask
setenv sfwf_filename               ${grid}_sfwf
setenv shf_filename                ${grid}_shf
setenv tidal_mixing_filename       ${grid}_tidal_mixing
setenv topography_filename         ${grid}_topography
setenv vert_grid_filename          ${grid}_vert_grid
setenv transport_contents_filename ${grid}_transport_contents
setenv tavg_contents_filename      ${grid}_tavg_contents
setenv history_contents_filename   ${grid}_history_contents
setenv movie_contents_filename     ${grid}_movie_contents

setenv DIN_LOC_ROOT  /fis/cgd/cseg/csm/inputdata

if ($grid == gx1v5) then
    setenv init_ts_filename        $DIN_LOC_ROOT/ocn/pop/gx1v5/ic/ts_PHC2_jan_ic_gx1v5_20061230.ieeer8
    setenv horiz_grid_filename     $DIN_LOC_ROOT/ocn/pop/gx1v5/grid/horiz_grid_20010402.ieeer8 
    setenv regionmask_filename     $DIN_LOC_ROOT/ocn/pop/gx1v5/grid/region_mask_20061229.ieeei4 
    setenv topography_filename     $DIN_LOC_ROOT/ocn/pop/gx1v5/grid/topography_20061229.ieeei4 
    setenv shf_filename            $DIN_LOC_ROOT/ocn/pop/gx1v5/forcing/shf_mm_all_85-88_20010308.ieeer8 
    setenv sfwf_filename           $DIN_LOC_ROOT/ocn/pop/gx1v5/forcing/sfwf_mm_PHC2_salx_flxio_20061230.ieeer8 
    setenv tidal_mixing_filename   $DIN_LOC_ROOT/ocn/pop/gx1v5/forcing/tidal_energy_gx1v5_20070102.ieeer8 
    setenv chl_filename            $DIN_LOC_ROOT/ocn/pop/gx1v5/forcing/chl_filled_gx1v5_20061230.ieeer8 
else if ($grid ==  gx3v5) then
    setenv init_ts_filename        $DIN_LOC_ROOT/ocn/pop/gx3v5/ic/ts_PHC2_jan_20030806.ieeer8
    setenv horiz_grid_filename     $DIN_LOC_ROOT/ocn/pop/gx3v5/grid/horiz_grid_20030806.ieeer8 
    setenv regionmask_filename     $DIN_LOC_ROOT/ocn/pop/gx3v5/grid/region_mask_20040220.ieeei4 
    setenv topography_filename     $DIN_LOC_ROOT/ocn/pop/gx3v5/grid/topography_20040323.ieeei4 
    setenv shf_filename            $DIN_LOC_ROOT/ocn/pop/gx3v5/forcing/shf_20031208.ieeer8 
    setenv sfwf_filename           $DIN_LOC_ROOT/ocn/pop/gx3v5/forcing/sfwf_20040517.ieeer8 
    setenv chl_filename            $DIN_LOC_ROOT/ocn/pop/gx3v5/forcing/chl_mm_SeaWiFs97-01_20031205.ieeer8 
else
    setenv init_ts_filename        unknown_init_ts
    setenv horiz_grid_filename     unknown_horiz_grid
    setenv regionmask_filename     unknown_regionmask
    setenv topography_filename     unknown_topography
    setenv shf_filename            unknown_shf
    setenv sfwf_filename           unknown_sfwf
    setenv chl_filename            unknown_chl
endif



./pop2_in_build.csh



end  # compset
end  # foreach OCN_GRID
