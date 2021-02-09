#!/bin/bash

#A Bash script to post-process (regriding, climatology generation and time-series extraction) to prepare E3SM model output to to be used in e3sm_diags.
#Note: For processing E3SM v2 output add '-P eam' or '-m eam' when calling ncclimo, also change 'cam' to 'eam' in this script

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified.sh
source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh

# Low-res Cori simulations

declare -a ids=("20180215.DECKv1b_H1.ne30_oEC.edison")
start='1980'
end='2014'

for caseid in "${ids[@]}"
do
echo ${caseid}

drc_in=/global/cfs/cdirs/e3smpub/E3SM_simulations/${caseid}/archive/atm/hist
drc_in_river=/global/cfs/cdirs/e3smpub/E3SM_simulations/${caseid}/archive/rof/hist
echo $drc_in

dir_out_climo=/global/cfs/cdirs/e3sm/zhang40/postprocessing_for_e3sm_diags/climo/${caseid}/${start}-${end}
dir_out_ts=/global/cfs/cdirs/e3sm/zhang40/postprocessing_for_e3sm_diags/monthly_ts/${caseid}/${start}-${end}
dir_out_diurnal_climo=/global/cfs/cdirs/e3sm/zhang40/postprocessing_for_e3sm_diags/diurnal_climo/${caseid}/${start}-${end}
map_file=/global/homes/z/zender/data/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc

echo "Generating Climatology files"
drc_rgr=${dir_out_climo}/rgr
drc_out=${dir_out_climo}/native
ncclimo --caseid=${caseid} --yr_srt=${start} --yr_end=${end} --drc_in=${drc_in} --drc_out=${drc_out} -O $drc_rgr --map=${map_file}

echo "Generating Diurnal Cycle Climo files"
drc_rgr=${dir_out_diurnal_climo}/rgr
drc_out=${dir_out_diurnal_climo}/native
echo ${drc_in}

cd ${drc_in};eval ls ${caseid}.cam.h4.*{${start}..${end}}*.nc | ncclimo --clm_md=hfc --caseid=${caseid}.cam.h4 -v PRECT --ypf=1 --yr_srt=${start} --yr_end=${end} --drc_out=${drc_out} -O $drc_rgr --map=${map_file}
##
echo "Generating per-variable monthly time-series."
drc_rgr=${dir_out_ts}/rgr
drc_out=${dir_out_ts}/native

echo "Variables for Streamflow"
cd ${drc_in_river};eval ls ${caseid}*mosart.h0.*{${start}..${end}}*.nc | ncclimo  --caseid=${caseid} --var_xtr=areatotal2 -v RIVER_DISCHARGE_OVER_LAND_LIQ --yr_srt=$start --yr_end=$end --drc_out=${drc_rgr}
#
echo "Variables for supporting diags using monthly time series as input(ENSO, QBO, etc.)"
cd ${drc_in};eval ls ${caseid}.cam.h0.*{${start}..${end}}*.nc | ncclimo  --caseid=${caseid} --var=U,CLDHGH,CLDLOW,CLDMED,CLDTOT,FLNS,FLUT,FSNS,FSNT,FSNTOA,LANDFRAC,LHFLX,LWCF,OCNFRAC,PRECC,PRECL,PSL,QFLX,SHFLX,SWCF,T,TAUX,TAUY,TREFHT,TS --yr_srt=$start --yr_end=$end --drc_out=${drc_out} -O $drc_rgr --map=${map_file}
done

exit
