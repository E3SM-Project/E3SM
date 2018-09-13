#!/bin/bash
# Script initiated by Jerry Potter and modified by Jill Zhang (zhang40@llnl.gov)

path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/MERRA2/'

original_data_path=$path'original_data/'
time_series_output_path=$path'time_series/'
climo_output_path=$path'climatology/'

mkdir $time_series_output_path
mkdir $climo_output_path

#start_yr=1979
start_yr=1980
end_yr=2015

cd $original_data_path
#Get varibles
#for file0 in *197901-*nc
for file0 in *198001-*nc
do 
    var=$(echo $file0 | cut -d'_' -f1)
    vars="$vars $var"

#    for file1 in ${var}_*nc
#    do
#        echo $file1
#        cdo splityear $file1 ${time_series_output_path}MERRA2_yearly_${var}_
#    done
#
#    for yr in $(eval echo "{$start_yr..$end_yr}"); do
#        echo "$yr"
#        yyyy=`printf "%04d" $yr`
#        for mth in {1..12}; do
#            mm=`printf "%02d" $mth`
#           ncks -O -F -d time,${mth} ${time_series_output_path}MERRA2_yearly_${var}_${yyyy}.nc ${time_series_output_path}MERRA2_monthly_${var}_${yyyy}${mm}.nc
#        done
#    done
#    ncrcat ${time_series_output_path}MERRA2_yearly_${var}_${yyyy}.nc ${time_series_output_path}MERRA2_${var}_${start_yr}_${end_yr}.nc 
#    ncclimo -a sdd --lnk_flg -c MERRA2_monthly_${var}_${start_yr}01.nc -s $start_yr -e $end_yr -i ${time_series_output_path} -o ${climo_output_path}
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
    ncks -A MERRA2_monthly_${i}_${j}_*nc MERRA2_monthly_vas_${j}_*nc
    done
done

for i in MERRA2_monthly_vas*.nc; do mv "$i" "${i/MERRA2_monthly_vas/MERRA2}" ; done
rm MERRA2_monthly_*.nc

