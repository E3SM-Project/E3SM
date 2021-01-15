include(ScreamUtils)

# Required libraries
set (NEED_LIBS rrtmgp scream_control scream_share yakl)

set (SRC rrtmgp_tests.cpp rrtmgp_test_utils.cpp)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
    set_source_files_properties(${SRC} PROPERTIES CXX_STANDARD 14)
    set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS "${YAKL_CXX_FLAGS}")
    CreateUnitTest(rrtmgp_tests "${SRC}" "${NEED_LIBS}" LABELS "rrtmgp")
    set_target_properties(rrtmgp_tests PROPERTIES CXX_STANDARD 14)
    target_include_directories(rrtmgp_tests PUBLIC ${SCREAM_BASE_DIR}/../../externals/YAKL)
endif()

# Copy RRTMGP absorption coefficient lookup tables to local data directory
# TODO: RRTMGP code and data should live in higher level "externals" directory
# not in the "eam" directory. We ultimately do not want dependencies on
# anything in EAM.
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/src/physics/rrtmgp/tests/rrtmgp-allsky.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
CONFIGURE_FILE(${SCREAM_BASE_DIR}/../eam/src/physics/rrtmgp/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)

set (SRC generate_baseline.cpp rrtmgp_test_utils.cpp)
CreateUnitTest(
    generate_baseline "${SRC}" "${NEED_LIBS}" LABELS "rrtmgp"
    EXCLUDE_MAIN_CPP
)
    #add_executable(
    #    generate_baseline "${SRC}" #"${NEED_LIBS}"
    #)
#target_link_libraries(generate_baseline ${NEED_LIBS})
#target_include_directories(generate_baseline PUBLIC ${SCREAM_BASE_DIR}/../../externals/YAKL)

add_custom_target(rrtmgp_allsky_baseline.nc
    COMMAND ${CMAKE_COMMAND} -E env $<TARGET_FILE:generate_baseline>)# ./generate_baseline)

add_dependencies(baseline rrtmgp_allsky_baseline.nc)
