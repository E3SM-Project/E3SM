#!/bin/bash

# Wrapper script for BGQ runjob command - used by the testing scripts
# This script is required because CMake ADD_TEST() does not provide a
# way to expand $COBALT_PARTNAME at runtime
args="$*"
runjob --block $COBALT_PARTNAME --verbose=OFF $@
