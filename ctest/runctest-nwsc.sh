#!/bin/sh
#==============================================================================
#
#  This script defines how to run CTest on the NCAR Wyoming Supercomputing 
#  Center (yellowstone/caldera/geyser) system.
#
#  This assumes the CTest model name (e.g., "Nightly") is passed to it when
#  run.
#
#==============================================================================

# Get the CTest model name
model=$1

# Run the "ctest" command through an interactive parallel session
DAV_CORES=4 execca ctest -S ctest/CTestScript-Test.cmake,${model} -V
