#!/bin/bash -e

## Set necesssary paths
export ALMUQ_PATH=${PWD}
export UQTK_PATH=${PWD}/UQTk_v2.2/src_cpp/
export UQTK_BIN=${UQTK_PATH}/bin

## If ncdump is on your path
export NCDUMP_PATH=

## On Titan, use 
if [[ $HOSTNAME == *"titan"* ]]; then
export NCDUMP_PATH=/opt/cray/netcdf/4.3.2/bin/
##and load matplotlib, too
module load python_matplotlib/1.2.1
fi
