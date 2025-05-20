#!/bin/bash

#export LD_PRELOAD=/usr/lib64/libpil4dfs.so
#printenv
echo "============== MPIEXEC for AURORA ==================="
echo "Launching job after pre-loading the DAOS Interception library..."
mpiexec -genv LD_PRELOAD=/usr/lib64/libpil4dfs.so $@
