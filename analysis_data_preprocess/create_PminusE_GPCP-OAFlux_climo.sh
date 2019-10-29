#!/bin/bash
# Script generated to process corev2_flux data into climatology and time_series files as input of e3sm_diags by Jill Zhang (zhang40@llnl.gov)

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/GPCP_OAFLux/'

time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
original_path=$path'original_data/'
tmp=$path'tmp/'

mkdir $climo_output_path
mkdir $time_series_output_path
mkdir $tmp

cd $time_series_output_path
start_yr=1979
end_yr=2013

cdo splityear PminusE_197901_201312.nc ${tmp}GPCP_OAFLux

for yr in $(eval echo "{$start_yr..$end_yr}"); do
    yyyy=`printf "%04d" $yr`
   

    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}GPCP_OAFLux${yyyy}.nc ${tmp}GPCP_OAFLux_${yyyy}${mm}.nc
        done
done
cd ${tmp}

ncclimo -a sdd --lnk_flg -c GPCP_OAFLux_${start_yr}01.nc -s $start_yr -e $end_yr
mv *climo.nc $climo_output_path









