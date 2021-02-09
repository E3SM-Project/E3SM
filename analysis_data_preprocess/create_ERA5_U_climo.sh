#!/bin/bash

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5_U/'

original_data_path=$path'original/ua/'

time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
tmp=$path'tmp/'

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp


start_yr=1979
end_yr=2018

dst_fl=/p/user_pub/e3sm/zhang40/destination_file/180_360.nc
drc_out=$tmp

for file0 in ${original_data_path}adaptor*nc
do
    echo $file0 | cut -d'-' -f9

    filename=${file0##*-} #retain the part after the last -
    echo $filename
    #switch latitude to N-to-S
    ncpdq -a time,-lat,lon $file0 ${tmp}ua_N-to-S_$filename
    #uncompress variable
    ncpdq -U ${tmp}ua_N-to-S_$filename ${tmp}ua_$filename
    #regrid from 0.25 deg to 1deg, convert double to single precison
    ncremap --d2f -d $dst_fl ${tmp}ua_$filename ${tmp}ua_180_360_$filename
    ncks --mk_rec_dmn time ${tmp}ua_180_360_$filename ${tmp}ua_180_360_rec_dmn_$filename
done
ncrcat ${tmp}ua_180_360_rec_dmn*nc ${time_series_output_path}ua_${start_yr}01_${end_yr}12.nc
ncrename -v u,ua ${time_series_output_path}ua_${start_yr}01_${end_yr}12.nc
