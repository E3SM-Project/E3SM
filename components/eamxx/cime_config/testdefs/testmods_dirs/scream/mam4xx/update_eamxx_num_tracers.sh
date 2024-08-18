#!/bin/sh

#------------------------------------------------------
# MAM4xx adds additional tracers to the simulation
# Increase number of tracers for MAM4xx simulations
#------------------------------------------------------

# Additional MAM4xx tracers (MAM4xx adds 31 tracers)
ADDITIONAL_MAM4xx_TRACERS=31

# Original CMAKE options in env_build.xml
orig_cmake_opt=`./xmlquery --value SCREAM_CMAKE_OPTIONS`

# Extract the number of tracers
orig_tracer_num=$(echo $orig_cmake_opt | grep -oP 'SCREAM_NUM_TRACERS \K[0-9]+')

# Update number of tracers
new_tracer_num=$((orig_tracer_num + ADDITIONAL_MAM4xx_TRACERS))

# Form the new CMake options string by replacing the original number with the new number
new_cmake_opt=$(echo $orig_cmake_opt | sed "s/SCREAM_NUM_TRACERS $orig_tracer_num/SCREAM_NUM_TRACERS $new_tracer_num/")

# Update cmake options string
`./xmlchange  SCREAM_CMAKE_OPTIONS="$new_cmake_opt"`
