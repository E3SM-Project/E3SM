#!/bin/bash
# Script initiated by Jerry Potter and modified by Jill Zhang (zhang40@llnl.gov)

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/ERA-Interim/'

original_data_path=$path'original_data/'
time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
tmp=$path'tmp/'

start_yr=1979
end_yr=2016

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp


cd $original_data_path
#Get varibles
#Note: in original data sets, 3d variables are stored one year of data per file

for file0 in *197901-*nc
do 
    var=$(echo $file0 | cut -d'_' -f1)
    vars="$vars $var"

    for file1 in ${var}_*nc
    do
        echo $file1
        cdo splityear $file1 ${tmp}ERA-Interim_yearly_${var}_
    done

    for yr in $(eval echo "{$start_yr..$end_yr}"); do
        echo "$yr"
        yyyy=`printf "%04d" $yr`
        for mth in {1..12}; do
            mm=`printf "%02d" $mth`
           ncks -O -F -d time,${mth} ${tmp}ERA-Interim_yearly_${var}_${yyyy}.nc ${tmp}ERA-Interim_monthly_${var}_${yyyy}${mm}.nc
        done
    done
    ncrcat ${tmp}ERA-Interim_yearly_${var}_*.nc ${time_series_output_path}${var}_${start_yr}01_${end_yr}12.nc 
    ncclimo -a sdd --lnk_flg -c ERA-Interim_monthly_${var}_${start_yr}01.nc -s $start_yr -e $end_yr -i ${tmp} -o ${climo_output_path}
#rm ${time_series_output_path}*yearly*nc
#rm ${time_series_output_path}*monthly*nc
done

#Combine all variables in one climo file then rename
#vars=hfls hfss pr prw psl rlds rlus rlut rlutcs rsds rsdt rsus rsut rsutcs tauu tauv uas vas
#vars=hfls hfss pr prw psl rlds rlus rlut rlutcs rsds rsdt rsus rsut rsutcs tauu tauv uas

cd ${climo_output_path}
echo $vars

declare -a sn=("ANN" "DJF" "MAM" "JJA" "SON" "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12")
for j in "${sn[@]}"
do
    #for i in "${var[@]}"
    for i in $vars
    do
    ncks -A ERA-Interim_monthly_${i}_${j}_*nc ERA-Interim_monthly_vas_${j}_*nc
    done
done

for i in ERA-Interim_monthly_vas*.nc; do mv "$i" "${i/ERA-Interim_monthly_vas/ERA-Interim}" ; done
rm ERA-Interim_monthly_*.nc

