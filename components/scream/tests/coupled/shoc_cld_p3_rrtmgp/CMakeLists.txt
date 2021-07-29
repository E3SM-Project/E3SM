INCLUDE (ScreamUtils)

## Add CUDA flags for YAKL
if (CUDA_BUILD)
    set(ARCH "CUDA")
    set(YAKL_CXX_FLAGS "-D__USE_CUDA__ --expt-extended-lambda --expt-relaxed-constexpr ${CUDA_FLAGS} ${YAKL_CXX_FLAGS}")
endif()

## Create test
# Required libraries
find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
SET (NEED_LIBS shoc cld_fraction p3 scream_rrtmgp rrtmgp ${NETCDF_C} scream_control scream_share physics_share yakl)

# Add source files
set (SRC shoc_cld_p3_rrtmgp.cpp ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests/rrtmgp_test_utils.cpp)

# Test atmosphere processes.
CreateUnitTest(shoc_cld_p3_rrtmgp "${SRC}" "${NEED_LIBS}" LABELS "shoc;cld;p3;rrtmgp;physics")
set_target_properties(shoc_cld_p3_rrtmgp PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
target_include_directories(shoc_cld_p3_rrtmgp PUBLIC
    ${SCREAM_BASE_DIR}/src/physics/rrtmgp
    ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/extensions/cloud_optics
)

## Copy yaml input file to run directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
# P3 data files
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
# RRTMGP data files
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/examples/all-sky/rrtmgp-allsky.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
# Basic config files
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
               ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/shoc_cld_p3_rrtmgp_init_ne2np4.nc
               ${CMAKE_CURRENT_BINARY_DIR}/shoc_cld_p3_rrtmgp_init_ne2np4.nc COPYONLY)
