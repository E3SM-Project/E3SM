INCLUDE (ScreamUtils)

SET (NEED_LIBS cld_fraction scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(cld_fraction_stand_alone "cld_fraction_stand_alone.cpp" "${NEED_LIBS}" LABELS "cld_fraction;physics")

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)
