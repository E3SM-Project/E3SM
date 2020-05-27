#!/bin/bash

set -e

conda create -y -n test --override-channels -c conda-forge -c e3sm -c defaults \
    python=3.7 "compass=*=nompi*"
