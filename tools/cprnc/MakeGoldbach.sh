#!/bin/sh
export PATH=/usr/local/pgi/linux86-64/bin:$PATH

#export NETCDF=/usr/local/netcdf-pgi
export NETCDF=/usr/local/netcdf_c-4.3.0_f-4.4-beta1-pgi-pgcc-pghf-13.7

gmake
mv cprnc cprnc.exe
cp cprnc.sh cprnc
