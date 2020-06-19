INCLUDE (ScreamUtils)

# Required libraries
SET (rrtmgpLib rrtmgp)
SET (NEED_LIBS ${rrtmgpLib} scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(rrtmgp_stand_alone "rrtmgp_stand_alone.cpp" "${NEED_LIBS}" LABELS "rrtmgp;physics")

# Copy RRTMGP absorption coefficient lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../cam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
