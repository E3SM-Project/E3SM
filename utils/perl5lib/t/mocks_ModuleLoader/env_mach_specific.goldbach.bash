#!/usr/bin/env bash -f 
#===============================================================================
# Automatically generated module settings for goldbach
#===============================================================================

.  /usr/share/Modules/init/bash
module purge  
if [ "$COMPILER" = "intel" ]
then
	module load compiler/intel/14.0.2
fi
if [ "$COMPILER" = "pgi" ]
then
	module load compiler/pgi/14.10
fi
if [ "$COMPILER" = "nag" ]
then
	module load compiler/nag/5.3.1-907
fi
if [ "$COMPILER" = "gnu" ]
then
	module load compiler/gnu/4.4.7
fi
export P4_GLOBMEMSIZE=500000000
export NETCDF_DIR=$NETCDF_PATH
ulimit -s unlimited
ulimit -c unlimited
