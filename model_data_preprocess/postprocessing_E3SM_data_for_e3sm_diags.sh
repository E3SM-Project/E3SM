#!/bin/bash

#A Bash script to post-process (regriding, climatology generation and time-series extraction) to prepare E3SM model output to to be used in e3sm_diags.

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh

#For typical EAM v2 ne gp2 grids.
caseid="20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
#res = 30/120 for ne30/ne120 grids
res=30               # res = 30/120 for ne30/ne120 grids
pg2=true             # Set false for v1 or true for v2 production simulations
atm_name="eam"       # Use "cam" for v1 or "eam" for v2 production simulations

start='0051'
end='0060'
input_path=/global/cfs/cdirs/e3smdata/zppy_complete_run_nersc_output/${caseid}
result_dir=/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags/${caseid}
map_file=/global/homes/z/zender/data/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc

###For typical v1 ne SE grids.
#caseid="20180215.DECKv1b_H1.ne30_oEC.edison"
#
##res = 30/120 for ne30/ne120 grids
#res=30               # res = 30/120 for ne30/ne120 grids
#pg2=false             # Set false for v1 or true for v2 production simulations
#atm_name="cam"       # Use "cam" for v1 or "eam" for v2 production simulations
#result_dir=/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v1_data_for_e3sm_diags/${caseid}
#input_path=/global/cfs/cdirs/e3smpub/E3SM_simulations/${caseid}
#map_file=/global/homes/z/zender/data/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc
#
#start='1980'
#end='2014'


#for caseid in "${ids[@]}"
#do
echo ${caseid}

drc_in=${input_path}/archive/atm/hist
drc_in_river=${input_path}/archive/rof/hist
echo $drc_in



dir_out_climo=$result_dir/climatology/


dir_out_ts=$result_dir/time-series/
dir_out_diurnal_climo=$result_dir/diurnal_climatology


echo "Generating Climatology files"
drc_rgr=${dir_out_climo}/rgr
drc_out=${dir_out_climo}/native
ncclimo -P ${atm_name} --caseid=${caseid} --yr_srt=${start} --yr_end=${end} --drc_in=${drc_in} --drc_out=${drc_out} -O $drc_rgr --map=${map_file}

echo "Generating Diurnal Cycle Climo files"
drc_rgr=${dir_out_diurnal_climo}/rgr
drc_out=${dir_out_diurnal_climo}/native
echo ${drc_in}

cd ${drc_in};eval ls ${caseid}.${atm_name}.h4.*{${start}..${end}}*.nc | ncclimo -P ${atm_name} --clm_md=hfc --caseid=${caseid}.${atm_name}.h4 -v PRECT --yr_srt=${start} --yr_end=${end} --drc_out=${drc_out} -O $drc_rgr --map=${map_file}
##
echo "Generating per-variable monthly time-series."
drc_rgr=${dir_out_ts}/rgr
drc_out=${dir_out_ts}/native

echo "Variables for Streamflow"
cd ${drc_in_river};eval ls ${caseid}*mosart.h0.*{${start}..${end}}*.nc | ncclimo --caseid=${caseid} --var_xtr=areatotal2 -v RIVER_DISCHARGE_OVER_LAND_LIQ --yr_srt=$start --yr_end=$end --drc_out=${drc_rgr} --split

echo "Variables for supporting diags using monthly time series as input(ENSO, QBO, mixed-phase partition etc.)"
cd ${drc_in};eval ls ${caseid}.${atm_name}.h0.*{${start}..${end}}*.nc | ncclimo -P ${atm_name} --caseid=${caseid} --var=U,CLDICE,CLDLIQ,CLDHGH,CLDLOW,CLDMED,CLDTOT,FLNS,FLUT,FSNS,FSNT,FSNTOA,LANDFRAC,LHFLX,LWCF,OCNFRAC,PRECC,PRECL,PSL,QFLX,SHFLX,SWCF,T,TAUX,TAUY,TREFHT,TS --yr_srt=$start --yr_end=$end --drc_out=${drc_out} -O $drc_rgr --map=${map_file} --split
#done

exit
