#!/bin/bash

singularity pull e3sm.sif \
    docker://ghcr.io/mahf708/e3sm-imgs:0.0.3

if [ $? -ne 0 ]; then
    exit -1
fi
