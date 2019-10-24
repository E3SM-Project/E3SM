#!/bin/bash
# Script generated to process corev2_flux data into climatology and time_series files as input of e3sm_diags by Jill Zhang (zhang40@llnl.gov)

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/COREv2_Flux/'

time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
original_path=$path'original_data/'
tmp=$path'tmp/'

mkdir $climo_output_path
mkdir $time_series_output_path
mkdir $tmp

cd $original_path
echo $path

for yr in {1979..2006}; do
    yyyy=`printf "%04d" $yr`
    ncks --mk_rec_dmn time ${original_path}${yyyy}.nc ${tmp}time_rec_dim_${yyyy}.nc

    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}time_rec_dim_${yyyy}.nc ${tmp}COREv2_Flux_${yyyy}${mm}.nc
        done
done

ncrcat ${tmp}COREv2_Flux_*nc ${time_series_output_path}COREv2_Flux_197901_200612.nc 

ncrename -v F_evap,evspsbl -v F_prec,pr -v F_roff,mrro -v Q_lat,hfls -v Q_sen,hfss -v Q_lwdn,rlds -v Q_lwup,rlus -v Q_swnet,rss -v taux,tauu -v tauy,tauv ${time_series_output_path}COREv2_Flux_197901_200612.nc

for var in evspsbl pr mrro hfls hfss rlds rlus rss tauu tauv
do
   echo $var
   ncks -v $var ${time_series_output_path}COREv2_Flux_197901_200612.nc ${time_series_output_path}${var}_197901_200612.nc
done
rm ${time_series_output_path}COREv2_Flux*nc







