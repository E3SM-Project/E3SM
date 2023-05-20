#!/bin/bash

# This script is to generate climatology and time-series based on ERA5 variables not included in obs4mip archive. pr and et are duplicated for cross-validation.
path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5/'

## Set up for one extension
#ext='ext_1'
## variables include:
## si10: 10 metre wind speed
## d2m: 2 metre dewpoint temperature
## sp: Surface pressure
#
#declare -a var=("sp" "d2m" "si10")
#filename='adaptor.mars.internal-1683916975.6853056-27160-9-445faae7-c61c-4ee6-8cb2-32d1c02b0839.nc'


ext='ext'
# variables include:
# t2m: 2 meter temp"
# cp: Convective precipitation
# e: Evaporation
# lsp: Large-scale precipitation
# ro: Runoff
# tp: Total precipitation
# vimd: Vertically integrated moisture divergence
declare -a var=("t2m" "cp" "e" "lsp" "ro" "tp" "vimd")
filename='adaptor.mars.internal-1649447170.9796565-18358-9-69c2693a-cb8f-49a8-8f09-36e7e03c7239.nc'
original_data_path=$path'original_'${ext}'/'
time_series_output_path=$path'time_series_'${ext}'_d2f/'
climo_output_path=$path'climatology_'${ext}'_d2f/'
tmp=$path'tmp_'${ext}'_d2f/'

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp

start_yr=1979
end_yr=2019

## reduce expver dimension (only needed for recent data), reference: https://confluence.ecmwf.int/display/CUSF/ERA5+CDS+requests+which+return+a+mixture+of+ERA5+and+ERA5T+data 
#cdo --reduce_dim -copy ${original_data_path}${filename} ${tmp}reduce_expver.nc

#switch latitude to N-to-S
ncpdq -a time,-lat,lon ${original_data_path}${filename} ${tmp}N-to-S.nc
#unpack variable then convert from double to float
ncpdq -O --unpack ${tmp}N-to-S.nc ${tmp}N-to-S_unpack.nc
ncpdq -O --map=dbl_flt ${tmp}N-to-S_unpack.nc ${tmp}N-to-S_unpack_d2f.nc
#fix attribute
ncatted --attribute=missing_value,,d,, --attribute=_FillValue,,d,, ${tmp}N-to-S_unpack_d2f.nc

cdo splityear ${tmp}N-to-S_unpack_d2f.nc ${tmp}ERA5_ext

for yr in $(eval echo "{$start_yr..$end_yr}"); do
    yyyy=`printf "%04d" $yr`
    echo $yyyy

    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}ERA5_ext${yyyy}.nc ${tmp}ERA5_ext_${yyyy}${mm}.nc
        done
done

cd ${tmp}

ncclimo -a sdd -c ${tmp}ERA5_ext_${start_yr}01.nc -s $start_yr -e $end_yr
mv *climo.nc $climo_output_path

ncrcat ${tmp}ERA5_ext_*nc ${time_series_output_path}ERA5_ext_${start_yr}01_${end_yr}12.nc

# time series of variables are splitted into one variable each file, ex:
# declare -a var=("sp" "d2m" "si10")
for i in "${var[@]}"
do 
    ncks -v ${i} ${time_series_output_path}ERA5_ext_${start_yr}01_${end_yr}12.nc ${time_series_output_path}${i}_${start_yr}01_${end_yr}12.nc
done
mv ${time_series_output_path}ERA5_ext_${start_yr}01_${end_yr}12.nc ${tmp}

# climatology are appended
declare -a sn=("ANN" "DJF" "MAM" "JJA" "SON" "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
for j in "${sn[@]}"
do
    ncks -A /p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5/climatology_ext_1_d2f/ERA5_ext_${j}_*nc /p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5/climatology_ext_d2f/ERA5_ext_${j}_*nc
done


exit
