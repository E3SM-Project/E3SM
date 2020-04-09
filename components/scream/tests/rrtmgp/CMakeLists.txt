INCLUDE (ScreamUtils)

#SET (p3Lib p3)
#SET (NEED_LIBS ${p3Lib} scream_control scream_share)
#SET (rrtmgpLib rrtmgp)
SET (NEED_LIBS rrtmgp scream_control scream_share)

# Test atmosphere processes
CreateUnitTest(rrtmgp_stand_alone "rrtmgp_stand_alone.cpp" "${NEED_LIBS}")

# Copy lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
