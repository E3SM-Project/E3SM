#!/bin/bash

#SBATCH  --job-name=e3sm_diags_{{ model }}
#SBATCH  --nodes=1
#SBATCH  --output=e3sm_diags_{{ model }}.o%j
#SBATCH  --exclusive
#SBATCH  --time=02:00:00

# Load environment
#source /export/golaz1/conda/etc/profile.d/conda.sh
#conda activate e3sm_diags_env_dev
source /p/user_pub/e3sm_unified/envs/load_latest_e3sm_unified_acme1.sh

# Make sure UVCDAT doesn't prompt us about anonymous logging
export UVCDAT_ANONYMOUS_LOG=False

# Run E3SM Diags
time python << EOF

import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

param.test_data_path = '{{ simulation }}'

param.short_test_name = '{{ institution }} {{ model }} {{ experiment }} ({{ realization }})'
param.test_timeseries_input = True
param.test_start_yr = '1985'
param.test_end_yr = '2014'

#param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
param.reference_data_path = '/p/user_pub/e3sm/diagnostics/observations/Atm/climatology_1985-2014'


param.results_dir = '/var/www/acme/acme-diags/zhang40/CMIP6_20240109_1985-2014/{{ model }}/{{ experiment }}/{{ realization }}'
param.multiprocessing = True
param.num_workers = 16
#param.output_format_subplot = ["pdf", "png"]
param.diff_title = '{{ model }} {{ experiment }} ({{ realization }}) vs Obs'

# Use below to run all core sets of diags:
#runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
# Use below to run only a subset:
runner.sets_to_run = ['lat_lon']#, 'zonal_mean_xy', 'zonal_mean_2d']
runner.run_diags([param])

EOF

echo All done...


