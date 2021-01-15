INCLUDE (ScreamUtils)

# Required libraries
set (NEED_LIBS rrtmgp scream_control scream_share yakl)

set (SRC rrtmgp_stand_alone.cpp ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests/rrtmgp_test_utils.cpp)

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

# Test atmosphere processes
if (NOT ${SCREAM_BASELINES_ONLY})
    CreateUnitTest(rrtmgp_stand_alone "${SRC}" "${NEED_LIBS}" LABELS "rrtmgp;physics")
    set_source_files_properties("${SRC}" PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
    target_include_directories(rrtmgp_stand_alone PUBLIC ${SCREAM_BASE_DIR}/../../externals/YAKL)
    target_include_directories(rrtmgp_stand_alone PUBLIC ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests)
endif()

# Copy yaml input file to run directory
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

# Copy RRTMGP absorption coefficient lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/rad_standalone_ne4_pg2.nc
               ${CMAKE_CURRENT_BINARY_DIR}/rad_standalone_ne4_pg2.nc COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests/rrtmgp-allsky.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
