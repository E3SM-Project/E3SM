INCLUDE (ScreamUtils)

SET (cldfracLib cldfrac)
SET (NEED_LIBS ${cldfracLib} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(cldfrac_stand_alone "cld_fraction_stand_alone.cpp" "${NEED_LIBS}" LABELS "cldfrac;physics")

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

CONFIGURE_FILE(cld_fraction_standalone_ne4_ne4.nc cld_fraction_standalone_ne4_ne4.nc COPYONLY)
CONFIGURE_FILE(cld_fraction_standalone_output.yaml cld_fraction_standalone_output.yaml COPYONLY)
CONFIGURE_FILE(cld_fraction_standalone_ne4_ne4.baseline.nc cld_fraction_standalone_ne4_ne4.baseline.nc COPYONLY)
