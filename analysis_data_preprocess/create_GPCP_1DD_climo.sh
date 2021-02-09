#! /usr/bin/env bash

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/GPCP_Daily_1DD/'

original_data_path=$path'original_data/'
time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
tmp=$path'tmp/'

#start_yr=1997
#end_yr=2016

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp


cd $original_data_path
for yyyy in {1997..2016}; do      # Loop over years
  for moy in {1..12}; do          # Loop over months
    mm=$( printf "%02d" ${moy} )  # Change to 2-digit format

    # Average specific month yyyy-mm
    ncra -O -d time,"${yyyy}-${mm}-01","${yyyy}-${mm}-31" \
         ${original_data_path}pr_day_GPCP-1-3_BE_gn_19961002-20170101.nc ${tmp}GPCP_1DD_${yyyy}${mm}.nc
  done
done
#ncrcat ${tmp}GPCP_1DD_*.nc ${time_series_output_path}pr_199701_201612.nc
ncclimo -a sdd --no_amwg_link -c GPCP_1DD_199701.nc -s 1997 -e 2016 -i ${tmp} -o ${climo_output_path}

exit

# Concatenate monthly files together
#ncrcat -O in_??????.nc out.nc
