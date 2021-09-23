#!/bin/bash

wget -t 3 -O e3sm.sif http://portal.nersc.gov/project/e3sm/lukasz/e3sm.sif || \
    wget -t 3 -O e3sm.sif https://dabdceba-6d04-11e5-ba46-22000b92c6ec.e.globus.org/containers/public/e3sm.sif
if [ $? -ne 0 ]; then
    exit -1
fi
