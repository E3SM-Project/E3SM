#!/bin/bash
source /usr/local/e3sm_unified/envs/load_latest_e3sm_unified_acme1.sh

# Low-res Cori simulations
#drc_in=/p/user_pub/PCMDIobs/PCMDIobs2/atmos/3hr/pr/TRMM-3B43v-7/gn/v20200707
drc_in=/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/TRMM/original

drc_out=/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/TRMM/climatology_diurnal_cycle
caseid=pr_3hr_TRMM-3B43v-7_BE_gn_v20200707

cd ${drc_in};ls ${caseid}*.nc | ncclimo --var=pr --clm_md=hfc --caseid=TRMM-3B43v-7_3hr --ypf=1 --yr_srt=1998 --yr_end=2013 --drc_out=${drc_out}
