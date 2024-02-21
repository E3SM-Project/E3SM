#!/bin/bash

singularity pull e3sm.sif \
    docker://ghcr.io/mahf708/e3sm-imgs@sha256:d1030a6f4e3a53f682859436a26b30a9477d69423829ae1d9c1b5ab4e255430d

if [ $? -ne 0 ]; then
    exit -1
fi
