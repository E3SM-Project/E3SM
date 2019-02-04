#!/bin/bash
# Script initiated by Jerry Potter and modified by Jill Zhang (zhang40@llnl.gov)

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/HadISST/'

original_data_path=$path'original_data/'
time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
tmp=$path'tmp/'


mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp

cd $original_data_path
echo $path

cdo splityear HadISST_sst_ice_masked.nc ${tmp}HadISST1

for yr in {1870..2017}; do
    yyyy=`printf "%04d" $yr`
#    cdo invertlat 'HadISST0'${yyyy}'.nc' 'HadISST1'${yyyy}'.nc'
    cdo setrtomiss,-1e20,-1000 ${tmp}'HadISST1'${yyyy}'.nc' ${tmp}'HadISST2'${yyyy}'.nc'
    rm ${tmp}'HadISST1'${yyyy}'.nc'

    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}HadISST2${yyyy}.nc ${tmp}HadISST_PI_${yyyy}${mm}.nc
        cp ${tmp}HadISST_PI_${yyyy}${mm}.nc ${tmp}HadISST_PD_${yyyy}${mm}.nc
        cp ${tmp}HadISST_PI_${yyyy}${mm}.nc ${tmp}HadISST_CL_${yyyy}${mm}.nc
        done
done

ncrcat ${tmp}HadISST2*nc ${time_series_output_path}sst_187001_201712.nc 

#add time_bnds to time attr
#ncatted -O -a bounds,time,a,c,time_bnds ${time_series_output_path}sst_187001_201712.nc

cd ${tmp}
ncclimo -a sdd --lnk_flg -c HadISST_PI_187001.nc -s 1870 -e 1900
ncclimo -a sdd --lnk_flg -c HadISST_PD_187001.nc -s 1999 -e 2016
ncclimo -a sdd --lnk_flg -c HadISST_CL_187001.nc -s 1982 -e 2011
mv *climo.nc $climo_output_path


#ncclimo -a sdd --lnk_flg -c HadISST_PI_187001.nc -s 1999 -e 2016
#ncclimo -a sdd --lnk_flg -c HadISST_PI_187001.nc -s 1870 -e 1900
#rm 'HadiSST'${yyyy}'.nc'
#for yr in {1870..2016}; do
#    yyyy=`printf "%04d" $yr`
#    cdo setrtomiss,-1e20,-1000 'HadiSSTi'${yyyy}'.nc' 'HadiSSTii'${yyyy}'.nc'
#    rm 'HadiSSTi'${yyyy}'.nc'
#    mv 'HadiSSTii'${yyyy}'.nc' 'HadiSST'${yyyy}'.nc'
#done




