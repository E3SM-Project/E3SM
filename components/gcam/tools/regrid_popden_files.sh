#!/bin/bash

# Create multiple population density files for E3SM-GCAM at different resolutions
# One input argument determines resolution (see below)
# These are all the same format as the original 0.5x0.5 files, but regridded to the desired resolution
#    This format is what is used for the population density ELM files for the fire model, which are specified by namelist

# this script starts with the 0.5x0.5 ssp popluation density netcdf files
# all 5 ssp files are processed in the same way, using ncremap, but with different output names

# five output resolutions currently available:
# 0.9x1.25 (aka f09)
# 1.9x2.5 (aka f19)
# 0.5x0.5 (aka r05)
# 0.25x0.25 (aka r025)  
# 0.125x0.125 (aka r0125)
# note that no grid remapping is needed for 0.5x0.5 output

# The desired output resolution is selected by a single argument: 
# f09 = 0.9x1.25
# f19 = 1.9x2.5
# r05 = 0.5x0.5
# r025 = 0.25x0.25
# r0125 = 0.125x0.125

if [ "$#" != 1 ]; then
   echo "Usage: $0 <output resolution>"
   echo "Currently supported resolutions are: f09, f19, r05, r025, and r0215"
   exit
fi

RES=$1

date

# needed modules
#module load intel
#module load nco

# this gets what is needed also
source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

proc_dir='/compyfs/inputdata/iac/giac/iac2gcam/'

# the original data are at 0.5x0.5 res with origin at -180,-90; and corners aligned with these limits 
# this is also how the nomask scrip file is defined for 0.5x0.5

# can add more resolutions here

if [ $RES == 'f09' ]; then

   # f09

   # grid mapping
   map_file="/compyfs/inputdata/lnd/clm2/mappingdata/maps/0.9x1.25/map_0.5x0.5_nomask_to_0.9x1.25_nomask_aave_da_c121019.nc"

   # output nc file names - make sure they match the map file out resolution
   ssp1_file_out=${proc_dir}'elmforc.ssp1_hdm_0.9x1.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp2_file_out=${proc_dir}'elmforc.ssp2_hdm_0.9x1.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp3_file_out=${proc_dir}'elmforc.ssp3_hdm_0.9x1.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp4_file_out=${proc_dir}'elmforc.ssp4_hdm_0.9x1.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp5_file_out=${proc_dir}'elmforc.ssp5_hdm_0.9x1.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'

elif [ $RES == 'f19' ]; then

   # f19

   # grid mapping
   map_file="/compyfs/inputdata/lnd/clm2/mappingdata/maps/1.9x2.5/map_0.5x0.5_nomask_to_1.9x2.5_nomask_aave_da_c120709.nc"

   # output nc file names - make sure they match the map file out resolution
   ssp1_file_out=${proc_dir}'elmforc.ssp1_hdm_1.9x2.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp2_file_out=${proc_dir}'elmforc.ssp2_hdm_1.9x2.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp3_file_out=${proc_dir}'elmforc.ssp3_hdm_1.9x2.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp4_file_out=${proc_dir}'elmforc.ssp4_hdm_1.9x2.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp5_file_out=${proc_dir}'elmforc.ssp5_hdm_1.9x2.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'

elif [ $RES == 'r0125' ]; then

   # r0125

   # grid mapping
   map_file="/compyfs/inputdata/lnd/clm2/mappingdata/maps/0.125x0.125/map_0.5x0.5_nomask_to_0.125x0.125_nomask_aave_da_c241205.nc"

   # output nc file names - make sure they match the map file out resolution
   ssp1_file_out=${proc_dir}'elmforc.ssp1_hdm_0.125x0.125_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp2_file_out=${proc_dir}'elmforc.ssp2_hdm_0.125x0.125_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp3_file_out=${proc_dir}'elmforc.ssp3_hdm_0.125x0.125_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp4_file_out=${proc_dir}'elmforc.ssp4_hdm_0.125x0.125_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp5_file_out=${proc_dir}'elmforc.ssp5_hdm_0.125x0.125_simyr1850-2101_c'$(date +%Y%m%d)'.nc'

elif [ $RES == 'r05' ]; then

   # r05

   # no grid mapping
   map_file=""

   # output nc file names - make sure they match the map file out resolution
   ssp1_file_out=${proc_dir}'elmforc.ssp1_hdm_0.5x0.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp2_file_out=${proc_dir}'elmforc.ssp2_hdm_0.5x0.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp3_file_out=${proc_dir}'elmforc.ssp3_hdm_0.5x0.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp4_file_out=${proc_dir}'elmforc.ssp4_hdm_0.5x0.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp5_file_out=${proc_dir}'elmforc.ssp5_hdm_0.5x0.5_simyr1850-2101_c'$(date +%Y%m%d)'.nc'

elif [ $RES == 'r025' ]; then

   # r025

   # grid mapping
   map_file="/compyfs/inputdata/lnd/clm2/mappingdata/maps/0.25x0.25/map_0.5x0.5_nomask_to_0.25x0.25_nomask_aave_da_c250313.nc"

   # output nc file names - make sure they match the map file out resolution
   ssp1_file_out=${proc_dir}'elmforc.ssp1_hdm_0.25x0.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp2_file_out=${proc_dir}'elmforc.ssp2_hdm_0.25x0.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp3_file_out=${proc_dir}'elmforc.ssp3_hdm_0.25x0.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp4_file_out=${proc_dir}'elmforc.ssp4_hdm_0.25x0.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'
   ssp5_file_out=${proc_dir}'elmforc.ssp5_hdm_0.25x0.25_simyr1850-2101_c'$(date +%Y%m%d)'.nc'

else
	echo "$RES is not supported"
        echo "f09, f19, r0125, r025, and r05 are the currently supported output resolutions"
        exit
fi

#####

# source files
ssp1_source_file="/compyfs/inputdata/lnd/clm2/firedata/elmforc.ssp1_hdm_0.5x0.5_simyr1850-2101_c20200624.nc"
ssp2_source_file="/compyfs/inputdata/lnd/clm2/firedata/elmforc.ssp2_hdm_0.5x0.5_simyr1850-2101_c20200623.nc"
ssp3_source_file="/compyfs/inputdata/lnd/clm2/firedata/elmforc.ssp3_hdm_0.5x0.5_simyr1850-2101_c20200624.nc"
ssp4_source_file="/compyfs/inputdata/lnd/clm2/firedata/elmforc.ssp4_hdm_0.5x0.5_simyr1850-2101_c20200624.nc"
ssp5_source_file="/compyfs/inputdata/lnd/clm2/firedata/elmforc.ssp5_hdm_0.5x0.5_simyr1850-2100_c190109.nc"


# remap the data to the desired grid
# not needed for 0.5x0.5
if [ $RES != 'r05' ]; then
   ncremap -m ${map_file} ${ssp1_source_file} ${ssp1_file_out}
   ncremap -m ${map_file} ${ssp2_source_file} ${ssp2_file_out}
   ncremap -m ${map_file} ${ssp3_source_file} ${ssp3_file_out}
   ncremap -m ${map_file} ${ssp4_source_file} ${ssp4_file_out}
   ncremap -m ${map_file} ${ssp5_source_file} ${ssp5_file_out}
elif [ $RES == 'r05' ]; then
   cp ${ssp1_source_file} ${ssp1_file_out}
   cp ${ssp2_source_file} ${ssp2_file_out}
   cp ${ssp3_source_file} ${ssp3_file_out}
   cp ${ssp4_source_file} ${ssp4_file_out}
   cp ${ssp5_source_file} ${ssp5_file_out}
fi

date

