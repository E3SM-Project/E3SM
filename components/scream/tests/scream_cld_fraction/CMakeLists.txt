INCLUDE (ScreamUtils)

SET (cld_fractionLib cld_fraction)
SET (NEED_LIBS ${cld_fractionLib} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(cld_fraction_stand_alone "cld_fraction_stand_alone.cpp" "${NEED_LIBS}" LABELS "cld_fraction;physics")

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

CONFIGURE_FILE(cld_fraction_standalone_ne4_ne4.nc cld_fraction_standalone_ne4_ne4.nc COPYONLY)
CONFIGURE_FILE(cld_fraction_standalone_output.yaml cld_fraction_standalone_output.yaml COPYONLY)
CONFIGURE_FILE(cld_fraction_standalone_ne4_ne4.baseline.nc cld_fraction_standalone_ne4_ne4.baseline.nc COPYONLY)
