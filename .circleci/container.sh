#!/bin/bash

wget -t 3  http://portal.nersc.gov/project/e3sm/lukasz/e3sm.sif
if [ $? -eq 0 ]; then
    exit 0;
fi

wget -t 3 https://dabdceba-6d04-11e5-ba46-22000b92c6ec.e.globus.org/containers/public/e3sm.sif
if [ $? -eq 0 ]; then
    exit 0;
fi

exit -1
