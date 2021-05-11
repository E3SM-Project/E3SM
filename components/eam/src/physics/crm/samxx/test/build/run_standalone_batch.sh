#!/bin/bash
# 
#BSUB -P cli115
#BSUB -W 2:00
#BSUB -nnodes 1
#BSUB -q batch
#BSUB -J crm_standalone
#BSUB -o crm_standalone.%J
#BSUB -e crm_standalone.%J
#BSUB -u hannah6@llnl.gov
#BSUB -N

# To run this batch script:
# bsub run_standalone_batch.sh

# Load the python environment
# ( conda create --name crm_test_env --channel conda-forge netcdf4 numpy )
source activate crm_test_env

./runtest.sh
