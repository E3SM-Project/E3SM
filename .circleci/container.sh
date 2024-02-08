#!/bin/bash

singularity pull e3sm.sif \
    docker://ghcr.io/mahf708/e3sm-imgs@sha256:bb0facbb01b6221931e8702d7ad1427bc31a4bf0551cca3ee4561706ada39519

if [ $? -ne 0 ]; then
    exit -1
fi
