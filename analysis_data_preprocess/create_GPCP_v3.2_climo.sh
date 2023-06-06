#!/bin/bash
# Script generated to process GPCP v3.2 into climatology and time_series files as input of e3sm_diags by Jill Zhang (zhang40@llnl.gov)
# GPCP v3.2 related URL: https://earthdata.nasa.gov/esds/competitive-programs/measures/long-term-gpcp-precipitation" 

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/GPCP_v3.2/'

time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
original_path=$path'original/'
tmp=$path'tmp/'

mkdir $climo_output_path
mkdir $time_series_output_path
mkdir $tmp

cd $original_path
echo $path
start_yr=1983
end_yr=2021
##start_yr=1985
##end_yr=2014

#
for yr in $(eval echo "{$start_yr..$end_yr}"); do
    yyyy=`printf "%04d" $yr`
    echo $yyyy
    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        #assign time as record dimension
        ncks --mk_rec_dmn time ${original_path}GPCPMON_L3_${yyyy}${mm}_V3.2.nc4 ${tmp}time_rec_${yyyy}${mm}.nc
        #rotate longitude from [-180,180) to [0, 360).
        ncks -O --msa_usr_rdr -d lon,0.25,179.75 -d lon,-179.75,-0.25 ${tmp}time_rec_${yyyy}${mm}.nc ${tmp}GPCP_v3.2_${yyyy}${mm}.nc
        ncap2 -O -s 'where(lon < 0) lon=lon+360' ${tmp}GPCP_v3.2_${yyyy}${mm}.nc ${tmp}GPCP_v3.2_${yyyy}${mm}.nc
    done
done

cd ${tmp};eval ls GPCP_v3.2_{${start_yr}..${end_yr}}*.nc | ncclimo -a sdd -c GPCP_v3.2 -s $start_yr -e $end_yr --drc_out=${climo_output_path} 
#cd ${tmp};eval ls GPCP_v3.2_{1985..2014}*.nc | ncclimo -a sdd -c GPCP_v3.2 -s 1985 -e 2014 --drc_out=${climo_output_path} 

ncrcat ${tmp}GPCP_v3.2_*nc ${time_series_output_path}GPCP_v3.2_${start_yr}01_${end_yr}12.nc

exit


