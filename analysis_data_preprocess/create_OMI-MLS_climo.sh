#! /usr/bin/env bash
#Tropospheric column ozone (TCO) data is downloaded from https://acd-ext.gsfc.nasa.gov/Data_services/cloud_slice/new_data.html
#stratospheric column ozone (SCO) data is provided by J. R. Ziemke based on https://doi.org/10.5194/acp-19-3257-2019
#Both data are originally in ASCII format and converted into netCDF files by Qi Tang.

path='/Users/zhang40/Documents/ACME_simulations/O3/'

original_data_path=$path'original_data/'
time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
tmp=$path'tmp/'

#start_yr=1997
#end_yr=2016

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp


#Add lon dimention to zonel mean
ncap2 -s 'SCO=O3strat' ${original_data_path}O3strat_ZMK.nc ${time_series_output_path}SCO_200410_201712.nc
#cp ${original_data_path}O3strat_ZMK.nc ${time_series_output_path}SCO_200410_201712.nc
cdo splityear ${time_series_output_path}SCO_200410_201712.nc ${tmp}sco
cp ${original_data_path}tropO3clmn.nc ${time_series_output_path}TCO_200410_201712.nc
cdo splityear ${original_data_path}tropO3clmn.nc ${tmp}tco


cd $original_data_path
for yr in {2005..2017}; do      # Loop over years
    yyyy=`printf "%04d" $yr`
    echo $yyyy

    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}sco${yyyy}.nc ${tmp}sco_${yyyy}${mm}.nc
        ncks -O -F -d time,${mth} ${tmp}tco${yyyy}.nc ${tmp}OMI-MLS_${yyyy}${mm}.nc
        # It turned out attaching latidude depended SCO to TCO flipped latitude of TCO
        # Add ncpdq command to flip SCO lat first
        ncpdq -O -h -a -lat ${tmp}sco_${yyyy}${mm}.nc ${tmp}sco_${yyyy}${mm}.nc
        ncks -A ${tmp}sco_${yyyy}${mm}.nc ${tmp}OMI-MLS_${yyyy}${mm}.nc
        done
done
ncclimo -a sdd -c OMI-MLS_200501.nc -s 2005 -e 2017 -i ${tmp} -o ${climo_output_path}


exit

# Concatenate monthly files together
#ncrcat -O in_??????.nc out.nc
