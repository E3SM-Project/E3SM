#!/bin/sh
#==============================================================================
#
#  This script defines how to run CTest on the NCAR Wyoming Supercomputing
#  Center systems (yellowstone/caldera/geyser).
#
#  This assumes the CTest model name (e.g., "Nightly") is passed to it when
#  run.
#
#==============================================================================

# Get the CTest script directory
scrdir=$1

# Get the CTest model name
model=$2

# Run the "ctest" command through an interactive parallel session
DAV_CORES=4 execca ctest -S ${scrdir}/CTestScript-Test.cmake,${model} -V
