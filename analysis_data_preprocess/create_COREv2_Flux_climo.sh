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
start_yr=1979
end_yr=2006

for yr in $(eval echo "{$start_yr..$end_yr}"); do
    yyyy=`printf "%04d" $yr`
    echo $yyyy
    ncks --mk_rec_dmn time ${original_path}${yyyy}.nc ${tmp}time_rec_dim_${yyyy}.nc
    #ncrename -v F_evap,evspsbl -v F_prec,pr -v F_roff,mrro -v Q_lat,hfls -v Q_sen,hfss -v Q_lwdn,rlds -v Q_lwup,rlus -v Q_swnet,rss -v taux,tauu -v tauy,tauv ${original_path}${yyyy}.nc
    #somehow Q_lwup can not get renamed, neglect for now.
    #ncrename -v F_evap,evspsbl -v F_prec,pr -v F_roff,mrro -v Q_lat,hfls -v Q_sen,hfss -v Q_lwdn,rlds -v Q_swnet,rss -v taux,tauu -v tauy,tauv ${tmp}time_rec_dim_${yyyy}.nc


    for mth in {1..12}; do
        mm=`printf "%02d" $mth`
        ncks -O -F -d time,${mth} ${tmp}time_rec_dim_${yyyy}.nc ${tmp}COREv2_Flux_${yyyy}${mm}.nc
        done
done
cd ${tmp}

ncclimo -a sdd --lnk_flg -c COREv2_Flux_${start_yr}01.nc -s $start_yr -e $end_yr
mv *climo.nc $climo_output_path

ncrcat ${tmp}COREv2_Flux_*nc ${time_series_output_path}COREv2_Flux_${start_yr}01_${end_yr}12.nc 

#ncrename -v F_evap,evspsbl -v F_prec,pr -v F_roff,mrro -v Q_lat,hfls -v Q_sen,hfss -v Q_lwdn,rlds -v Q_lwup,rlus -v Q_swnet,rss -v taux,tauu -v tauy,tauv ${time_series_output_path}COREv2_Flux_${start_yr}01_${end_yr}12.nc

#for var in evspsbl pr mrro hfls hfss rlds rlus rss tauu tauv
for var in F_evap F_prec F_roff Q_lat Q_sen Q_lwdn Q_lwup Q_swnet taux tauy
do
   echo $var
   ncks -v $var ${time_series_output_path}COREv2_Flux_${start_yr}01_${end_yr}12.nc ${time_series_output_path}${var}_${start_yr}01_${end_yr}12.nc
done
rm ${time_series_output_path}COREv2_Flux*nc







