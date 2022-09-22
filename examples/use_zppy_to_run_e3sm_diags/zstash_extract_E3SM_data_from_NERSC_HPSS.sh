#!/usr/bin/env bash

# This script is an example to retrieve E3SM simulation data archived at NERSC HPSS to cfs file system for running E3SM Diags
# Note for extracting large data size, make sure to enter the data transfer node and to use "screen" before running this script
# For example, run following in the terminal. 
# $ ssh dtn01.nersc.gov
# $ screen

# Active e3sm_unified to use zstash
source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh

# Make a new folder to recieve extracted datasets.
mkdir -p /global/cfs/cdirs/e3smpub/E3SM_simulations/v2.LR.historical_0101
cd /global/cfs/cdirs/e3smpub/E3SM_simulations/v2.LR.historical_0101

# Note, the mapping between h0 to output frequency is defined in simulation config. Example below is valid for v2 LR Water Cycle Simulations
# Note, wildcards support is available in zstash extract (ig. to specify years to retrieve), refer to zstash usage here: https://e3sm-project.github.io/zstash/

# retrive h0(monthly) atmospheric data, size = 821 G
zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/v2.LR.historical_0101 "archive/atm/hist/*h0*"
# retrive h2(6 hourly) atmospheric data, size = 195 G [Optional] only needed for Tropical Cycle Analysis
zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/v2.LR.historical_0101 "archive/atm/hist/*h2*"
# retrive h4(3 hourly) atmospheric data, size = 195 G [Optional] only needed for Diurnal Cycle of Precipitation
zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/v2.LR.historical_0101 "archive/atm/hist/*h4*"
# retrive h0 (monthly) land data, size = 59 G [Optional]
zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/v2.LR.historical_0101 "archive/lnd/hist/*h0*"
# retrive h0 (monthly) river data, size = 59 G [Optional]
zstash extract --hpss=/home/projects/e3sm/www/WaterCycle/E3SMv2/LR/v2.LR.historical_0101 "archive/rof/hist/*h0*"

exit



