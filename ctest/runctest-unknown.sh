#!/bin/sh
#==============================================================================
#
#  This script defines how to run CTest on the default ("unknown") machine.
#
#  This assumes the CTest model name (e.g., "Nightly") is passed to it when
#  run.
#
#==============================================================================

# Get the dashboard model name
model=$1

# Run the "ctest" command in another process
ctest -S ctest/CTestScriptAppendix.cmake,${model} -I 5,8
