#!/bin/bash
# this script documents how to set up the landgen_env conda environment and install the landgen package in editable mode

# first load the conda module
module load conda

# first time only!
# create a the landgen_env conda environment from the .yml file
#conda env create -f landgen_env.yml

# to update the enviroment
conda env update -n landgen_env -f landgen_env.yml

# activate the landgen_env conda environment
conda activate landgen_env

# install a new python package
#conda install <package_name>