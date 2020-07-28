#!/bin/bash
# based on https://stackoverflow.com/a/49516361/7728169

set -e

source $HOME/miniconda/etc/profile.d/conda.sh
conda activate test

cd testing_and_setup/compass || exit 1
./list_testcases.py -h
./setup_testcase.py -h
./clean_testcase.py -h
./manage_regression_suite.py -h
cd ../..
