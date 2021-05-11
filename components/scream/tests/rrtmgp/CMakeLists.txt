INCLUDE (ScreamUtils)

# Add CUDA flags for YAKL
if (CUDA_BUILD)
    set(ARCH "CUDA")
    set(YAKL_CXX_FLAGS "-D__USE_CUDA__ --expt-extended-lambda --expt-relaxed-constexpr ${CUDA_FLAGS} ${YAKL_CXX_FLAGS}")
endif()
# Test atmosphere processes
if (NOT ${SCREAM_BASELINES_ONLY})

    # Required libraries
    find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
    set (NEED_LIBS scream_rrtmgp rrtmgp ${NETCDF_C} scream_control scream_share physics_share yakl)

    # Add source files
    set (SRC rrtmgp_stand_alone.cpp ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests/rrtmgp_test_utils.cpp)

    CreateUnitTest(
        rrtmgp_stand_alone "${SRC}" "${NEED_LIBS}" LABELS "rrtmgp;physics"
        EXE_ARGS "--ekat-test-params rrtmgp_inputfile=${CMAKE_CURRENT_BINARY_DIR}/data/rrtmgp-allsky.nc,rrtmgp_baseline=${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
    )
    message(STATUS "rrtmgp_stand_alone YAKL_CXX_FLAGS: ${YAKL_CXX_FLAGS}")
    set_target_properties(rrtmgp_stand_alone PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
    target_include_directories(rrtmgp_stand_alone PUBLIC
          ${SCREAM_BASE_DIR}/src/physics/rrtmgp
          ${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels
          ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/extensions/cloud_optics
    )

    # Copy yaml input file to run directory
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml
                   ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)

    # Copy RRTMGP initial condition to local data directory
    FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/rad_standalone_ne4_pg2.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/rad_standalone_ne4_pg2.nc COPYONLY)

    # Copy RRTMGP absorption coefficient lookup tables to local data directory
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
    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/rad_standalone_ne4_pg2.nc
                   ${CMAKE_CURRENT_BINARY_DIR}/rad_standalone_ne4_pg2.nc COPYONLY)
endif()
