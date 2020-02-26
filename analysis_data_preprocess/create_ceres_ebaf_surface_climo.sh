#!/bin/bash
## Script initiated by Jerry Potter and modified by Jill Zhang (zhang40@llnl.gov)
#directories yearly_data, monthly_data, and climo must be created. The original data is in the working directory.
#
path='/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/CERES-EBAF/'
cd $path


declare -a arr=("4.0" "2.8")
declare -a arr=("4.1")
for i in "${arr[@]}"
    do
        original_data_path=$path"$i"_surface/original_data/
        time_series_output_path=$path"$i"_surface/time_series/
        climo_data_output_path=$path"$i"_surface/climatology/
        mkdir $time_series_output_path
        mkdir $climo_data_output_path
    
        echo $path"$i"_surface/original_data/CERES_EBAF-Surface_Ed*.nc
        if [ $i == '4.1' ]
        then
            start_yr=2001
            end_yr=2018
            cdo chname,sfc_sw_down_all_mon,rsds,sfc_sw_down_clr_t_mon,rsdscs,sfc_sw_up_all_mon,rsus,sfc_sw_up_clr_t_mon,rsuscs,sfc_lw_down_all_mon,rlds,sfc_lw_down_clr_t_mon,rldscs,sfc_lw_up_all_mon,rlus,sfc_lw_up_clr_t_mon,rluscs ${original_data_path}CERES_EBAF-Surface_Ed*.nc ${original_data_path}ceres.nc
        else
            start_yr=2001
            end_yr=2015
            cdo chname,sfc_sw_down_all_mon,rsds,sfc_sw_down_clr_mon,rsdscs,sfc_sw_up_all_mon,rsus,sfc_sw_up_clr_mon,rsuscs,sfc_lw_down_all_mon,rlds,sfc_lw_down_clr_mon,rldscs,sfc_lw_up_all_mon,rlus,sfc_lw_up_clr_mon,rluscs ${original_data_path}CERES_EBAF-Surface_Ed*.nc ${original_data_path}ceres.nc
        fi

        echo 'Ed'"$i"
        cdo setmissval,1e20 ${original_data_path}ceres.nc ${original_data_path}ceres_miss.nc
        cdo splityear ${original_data_path}ceres_miss.nc ${time_series_output_path}ceres_ebaf_surface_time_series
        for yr in $(eval echo "{$start_yr..$end_yr}"); do
            echo "$yr"
            yyyy=`printf "%04d" $yr`
            for mth in {1..12}; do
                mm=`printf "%02d" $mth`
                echo ${mth}
               ncks -O -F -d time,${mth} ${time_series_output_path}ceres_ebaf_surface_time_series${yyyy}.nc ${time_series_output_path}ceres_ebaf_surface_${yyyy}${mm}.nc
            done
        done
        rm ${time_series_output_path}ceres_ebaf_surface_time_series*.nc
        ncrcat ${time_series_output_path}ceres_ebaf_surface_2*.nc ${time_series_output_path}ceres_ebaf_surface_v${i}_${start_yr}01_${end_yr}12.nc
        ncclimo -a sdd --lnk_flg -c ceres_ebaf_surface_${start_yr}01.nc -s $start_yr -e $end_yr -i ${time_series_output_path} -o ${climo_data_output_path}
        cd ${climo_data_output_path}
        for f in ceres_ebaf_surface_*.nc; do mv -v "$f" "${f/surface/surface_v$i}"; done;

#        rm ${time_series_output_path}ceres_ebaf_surface_v${i}_*nc
#        rm ${time_series_output_path}ceres_ebaf_surface_time_series*nc
   
done
