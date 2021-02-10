include(ScreamUtils)

# Required libraries
find_library(NETCDF_C netcdf HINTS ${NetCDF_C_PATHS}/lib)
set (NEED_LIBS rrtmgp ${NETCDF_C} scream_control scream_share yakl)

set (SRC rrtmgp_tests.cpp rrtmgp_test_utils.cpp)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
    set_source_files_properties(${SRC} PROPERTIES CXX_STANDARD 14)
    set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
    CreateUnitTest(
        rrtmgp_tests "${SRC}" "${NEED_LIBS}" LABELS "rrtmgp;physics"
        EXE_ARGS "${CMAKE_CURRENT_BINARY_DIR}/data/rrtmgp-allsky.nc ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
        EXCLUDE_MAIN_CPP
    )
    set_target_properties(rrtmgp_tests PROPERTIES CXX_STANDARD 14)
    target_include_directories(rrtmgp_tests PUBLIC ${SCREAM_BASE_DIR}/../../externals/YAKL)
endif()

# Build baseline code
set (GEN_BASELINE_SRC generate_baseline.cpp rrtmgp_test_utils.cpp)
add_executable(generate_baseline "${GEN_BASELINE_SRC}")
target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../../externals/YAKL)
target_include_directories(generate_baseline PUBLIC ${SCREAM_INCLUDE_DIRS})
target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp)
target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte)
target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rte/kernels)
target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp)
target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/rrtmgp/kernels)
target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/cpp/extensions/cloud_optics)
target_include_directories(generate_baseline PUBLIC tests)
set_target_properties(generate_baseline PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
set_target_properties(generate_baseline PROPERTIES CXX_STANDARD 14)
target_link_libraries(generate_baseline rrtmgp ${NETCDF_C} yakl)

# Copy RRTMGP absorption coefficient lookup tables to local data directory
# TODO: RRTMGP code and data should live in higher level "externals" directory
# not in the "eam" directory. We ultimately do not want dependencies on
# anything in EAM.
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
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

add_custom_target(rrtmgp_allsky_baseline.nc
    COMMAND ${CMAKE_COMMAND} -E env $<TARGET_FILE:generate_baseline>
    ${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/examples/all-sky/rrtmgp-allsky.nc ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc
)

add_dependencies(baseline rrtmgp_allsky_baseline.nc)
