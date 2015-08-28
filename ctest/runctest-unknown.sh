#!/bin/sh
#==============================================================================
#
#  This script defines how to run CTest on the default ("unknown") machine.
#
#  This assumes the CTest model name (e.g., "Nightly") is passed to it when
#  run.
#
#==============================================================================

# Run the "ctest" command in another process
ctest -I 5,6,7,8
