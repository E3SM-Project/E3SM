#!/bin/bash

# This script is to generate climatology and time-series based on ERA5 variables not included in obs4mip archive. pr and et are duplicated for cross-validation.
path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA5/'

original_data_path=$path'original_ext/'

time_series_output_path=$path'time_series_ext/'
climo_output_path=$path'climatology_ext/'
tmp=$path'tmp_ext/'

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp


start_yr=1979
end_yr=2019

# reduce expver dimension (only needed for recent data), reference: https://confluence.ecmwf.int/display/CUSF/ERA5+CDS+requests+which+return+a+mixture+of+ERA5+and+ERA5T+data 
#cdo --reduce_dim -copy ${original_data_path}adaptor.mars.internal-1649447170.9796565-18358-9-69c2693a-cb8f-49a8-8f09-36e7e03c7239.nc ${tmp}reduce_expver.nc

##switch latitude to N-to-S
#ncpdq -a time,-lat,lon ${original_data_path}adaptor.mars.internal-1649447170.9796565-18358-9-69c2693a-cb8f-49a8-8f09-36e7e03c7239.nc ${tmp}N-to-S.nc
##uncompress variable
#ncpdq -U ${tmp}N-to-S.nc ${tmp}N-to-S_uncompress.nc
##ncpdq -M dbl_flt
#
#cdo splityear ${tmp}N-to-S_uncompress.nc ${tmp}ERA5_ext

for yr in $(eval echo "{$start_yr..$end_yr}"); do
    yyyy=`printf "%04d" $yr`
    echo $yyyy
#    ncks --mk_rec_dmn time ${original_path}${yyyy}.nc ${tmp}time_rec_dim_${yyyy}.nc
    #ncrename -v F_evap,evspsbl -v F_prec,pr -v F_roff,mrro -v Q_lat,hfls -v Q_sen,hfss -v Q_lwdn,rlds -v Q_lwup,rlus -v Q_swnet,rss -v taux,tauu -v tauy,tauv ${original_path}${yyyy}.nc
    #somehow Q_lwup can not get renamed, neglect for now.
    #ncrename -v F_evap,evspsbl -v F_prec,pr -v F_roff,mrro -v Q_lat,hfls -v Q_sen,hfss -v Q_lwdn,rlds -v Q_swnet,rss -v taux,tauu -v tauy,tauv ${tmp}time_rec_dim_${yyyy}.nc

    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}ERA5_ext${yyyy}.nc ${tmp}ERA5_ext_${yyyy}${mm}.nc
        done
done

cd ${tmp}

ncclimo -a sdd -c ${tmp}ERA5_ext_${start_yr}01.nc -s $start_yr -e $end_yr
mv *climo.nc $climo_output_path

ncrcat ${tmp}ERA5_ext_*nc ${time_series_output_path}ERA5_ext_${start_yr}01_${end_yr}12.nc



#drc_out=$tmp
#
#for file0 in ${original_data_path}adaptor*nc
#do
#    echo $file0 | cut -d'-' -f9
#
#    #switch latitude to N-to-S
#    ncpdq -a time,-lat,lon $file0 ${tmp}ua_N-to-S_$filename
#    #uncompress variable
#    ncpdq -U --d2f ${tmp}ua_N-to-S_$filename ${tmp}ua_$filename
#    #regrid from 0.25 deg to 1deg, convert double to single precison
#    ncremap --d2f -d $dst_fl ${tmp}ua_$filename ${tmp}ua_180_360_$filename
#    ncks --mk_rec_dmn time ${tmp}ua_180_360_$filename ${tmp}ua_180_360_rec_dmn_$filename
#done
#ncrcat ${tmp}ua_180_360_rec_dmn*nc ${time_series_output_path}ua_${start_yr}01_${end_yr}12.nc
#ncrename -v u,ua ${time_series_output_path}ua_${start_yr}01_${end_yr}12.nc
