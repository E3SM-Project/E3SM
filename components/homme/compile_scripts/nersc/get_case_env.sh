#!/bin/bash
# Gets the case environment for standalone HOMME based on the `get_case_env`
# script of cime/CIME/Tools.
# Author: Jason Torchinsky

# Declare directory paths
# E3SM source code directory
e3sm=$(pwd)/../../../..

# Load the necessary modules
echo "-- Loading necessary modules..."
export e3sm
eval `$e3sm/cime/CIME/Tools/get_case_env`
echo "-- Loaded necessary modules!"

# Unset declared variables
unset e3sm

