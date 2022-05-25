#!/bin/bash
# Script generated to process aod550 MACv2 data into climatology and time_series files as input of e3sm_diags by Jill Zhang (zhang40@llnl.gov)
# MACv2 reference paper: https://www.tandfonline.com/doi/full/10.1080/16000889.2019.1623639

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/MACv2/'

time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
original_path=$path'original/'
tmp=$path'tmp/'

mkdir $climo_output_path
mkdir $time_series_output_path
mkdir $tmp

cd $original_path
echo $path
start_yr=2005
end_yr=2005

ncks -v aod ${original_path}gt_t_00550nm.nc ${tmp}aod550.nc
ncrename -v aod,AOD_550 ${tmp}aod550.nc ${tmp}AOD_550_2005.nc
ncatted -O -a long_name,AOD_550,m,c,"AOD at 550nm" ${tmp}AOD_550_2005.nc
ncatted -O -h -a contact,global,o,c,"Data title: MACv2 aerosol climatology\nData reference: Kinne, Stefan. The MACv2 aerosol climatology. Tellus B: Chemical and Physical Meteorology 71.1 (2019): 1-21. \nData creator: Chengzhu (Jill) Zhang, LLNL\nData contact:E3SM-DATA-SUPPORT@LLNL.GOV\nCreation date:07 April 2022\nData script: create_AOD550_MACv2_climo.sh\n" ${tmp}AOD_550_2005.nc 

for mth in {1..12}; do
    mm=`printf "%02d" $mth`
    ncks -O -F -d time,${mth} ${tmp}AOD_550_2005.nc ${tmp}MACv2_2005${mm}.nc
    done


cd ${tmp}

ncclimo -a sdd -c ${tmp}MACv2_${start_yr}01.nc -s $start_yr -e $end_yr
mv *climo.nc $climo_output_path

ncrcat ${tmp}MACv2_*nc ${time_series_output_path}MACv2_${start_yr}01_${end_yr}12.nc
rm -r $tmp


