#!/bin/sh
export PATH=/usr/local/pgi/linux86-64/bin:$PATH

export NETCDF=/usr/local/netcdf-pgi

gmake
mv cprnc cprnc.exe
cp cprnc.sh cprnc
