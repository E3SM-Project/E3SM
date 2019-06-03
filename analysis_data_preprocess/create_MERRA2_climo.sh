#!/bin/bash
# Script initiated by Jerry Potter and modified by Jill Zhang (zhang40@llnl.gov)

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/MERRA2/'

original_data_path=$path'original_data/'
time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'
tmp=$path'tmp/'

#start_yr=1979
start_yr=1980
end_yr=2016

mkdir $time_series_output_path
mkdir $climo_output_path
mkdir $tmp

cd $original_data_path
#Get varibles
#Note: in original data sets, 3d variables are stored one year of data per file
#for file0 in *197901-*nc
for file0 in *198001-*nc
do 
    var=$(echo $file0 | cut -d'_' -f1)
    echo "$var"
    vars="$vars $var"

    for file1 in ${var}_*nc
    do
        echo $file1
        cdo splityear $file1 ${tmp}MERRA2_yearly_${var}_
    done

    for yr in $(eval echo "{$start_yr..$end_yr}"); do
        echo "$yr"
        yyyy=`printf "%04d" $yr`
        for mth in {1..12}; do
            mm=`printf "%02d" $mth`
           ncks -O -F -d time,${mth} ${tmp}MERRA2_yearly_${var}_${yyyy}.nc ${tmp}MERRA2_monthly_${var}_${yyyy}${mm}.nc
        done
    done
    ncrcat ${tmp}MERRA2_yearly_${var}_*.nc ${time_series_output_path}${var}_${start_yr}01_${end_yr}12.nc 
    ncclimo -a sdd --lnk_flg -c MERRA2_monthly_${var}_${start_yr}01.nc -s $start_yr -e $end_yr -i ${tmp} -o ${climo_output_path}
#rm ${tmp}*yearly*nc
#rm ${tmp}*monthly*nc
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
    ncks -A MERRA2_monthly_${i}_${j}_*nc MERRA2_monthly_vas_${j}_*nc
    done
done

for i in MERRA2_monthly_vas*.nc; do mv "$i" "${i/MERRA2_monthly_vas/MERRA2}" ; done
rm MERRA2_monthly_*.nc

