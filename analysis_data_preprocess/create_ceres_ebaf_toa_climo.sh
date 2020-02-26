#!/bin/bash
## Script initiated by Jerry Potter and modified by Jill Zhang (zhang40@llnl.gov)
#directories yearly_data, monthly_data, and climo must be created. The original data is in the working directory.
#
path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/CERES-EBAF/'
cd $path


#declare -a arr=("4.0" "2.8")
declare -a arr=("4.1" "4.0")
for i in "${arr[@]}"
    do
        original_data_path=$path"$i"_toa/original_data/
        time_series_output_path=$path"$i"_toa/time_series/
        climo_data_output_path=$path"$i"_toa/climatology/
        mkdir $time_series_output_path
        mkdir $climo_data_output_path
    
        echo $path"$i"_toa/original_data/CERES_EBAF-TOA_Ed*.nc
        if [ $i == '4.1' ]
        then
            start_yr=2001
            end_yr=2018
            cdo chname,toa_sw_all_mon,rsut,toa_lw_all_mon,rlut,toa_sw_clr_t_mon,rsutcs,toa_lw_clr_t_mon,rlutcs,solar_mon,rsdt ${original_data_path}CERES_EBAF-TOA_Ed*.nc ${original_data_path}ceres.nc
        else
            start_yr=2001
            end_yr=2015
            cdo chname,toa_sw_all_mon,rsut,toa_lw_all_mon,rlut,toa_sw_clr_mon,rsutcs,toa_lw_clr_mon,rlutcs,solar_mon,rsdt ${original_data_path}CERES_EBAF-TOA_Ed*.nc ${original_data_path}ceres.nc
        fi

        echo 'Ed'"$i"
        cdo setmissval,1e20 ${original_data_path}ceres.nc ${original_data_path}ceres_miss.nc
        cdo splityear ${original_data_path}ceres_miss.nc ${time_series_output_path}ceres_ebaf_toa_time_series
        for yr in $(eval echo "{$start_yr..$end_yr}"); do
            echo "$yr"
            yyyy=`printf "%04d" $yr`
            for mth in {1..12}; do
                mm=`printf "%02d" $mth`
                echo ${mth}
               ncks -O -F -d time,${mth} ${time_series_output_path}ceres_ebaf_toa_time_series${yyyy}.nc ${time_series_output_path}ceres_ebaf_toa_${yyyy}${mm}.nc
            done
        done
        rm ${time_series_output_path}ceres_ebaf_toa_time_series*.nc
        ncrcat ${time_series_output_path}ceres_ebaf_toa_2*.nc ${time_series_output_path}ceres_ebaf_toa_v${i}_${start_yr}01_${end_yr}12.nc
        ncclimo -a sdd --lnk_flg -c ceres_ebaf_toa_${start_yr}01.nc -s $start_yr -e $end_yr -i ${time_series_output_path} -o ${climo_data_output_path}
        cd ${climo_data_output_path}
        for f in ceres_ebaf_toa_*.nc; do mv -v "$f" "${f/toa/toa_v$i}"; done;

#        rm ${time_series_output_path}ceres_ebaf_toa_v${i}_*nc
#        rm ${time_series_output_path}ceres_ebaf_toa_time_series*nc
   
done
