include(ScreamUtils)

# Required libraries
set (NEED_LIBS scream_rrtmgp)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
    CreateUnitTest(
        rrtmgp_tests rrtmgp_tests.cpp "${NEED_LIBS}" LABELS "rrtmgp;physics"
        EXE_ARGS "-i ${CMAKE_CURRENT_BINARY_DIR}/data/rrtmgp-allsky.nc -b ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
        EXCLUDE_MAIN_CPP
    )
    CreateUnitTest(
        rrtmgp_unit_tests rrtmgp_unit_tests.cpp "${NEED_LIBS}" LABELS "rrtmgp;physics"
    )
endif()

# Build baseline code
add_executable(generate_baseline generate_baseline.cpp)
target_link_libraries(generate_baseline PUBLIC ${NEED_LIBS})

# Copy RRTMGP absorption coefficient lookup tables to local data directory
# TODO: RRTMGP code and data should live in higher level "externals" directory
# not in the "eam" directory. We ultimately do not want dependencies on
# anything in EAM.
file (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)

configure_file(${EAM_RRTMGP_DIR}/external/rrtmgp/data/rrtmgp-data-sw-g224-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${EAM_RRTMGP_DIR}/external/rrtmgp/data/rrtmgp-data-lw-g256-2018-12-04.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${EAM_RRTMGP_DIR}/external/examples/all-sky/rrtmgp-allsky.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${EAM_RRTMGP_DIR}/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-sw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
configure_file(${EAM_RRTMGP_DIR}/external/extensions/cloud_optics/rrtmgp-cloud-optics-coeffs-lw.nc
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)

add_custom_target(rrtmgp_allsky_baseline.nc
    COMMAND ${CMAKE_COMMAND} -E env $<TARGET_FILE:generate_baseline>
            ${EAM_RRTMGP_DIR}/external/examples/all-sky/rrtmgp-allsky.nc
            ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc
)

add_dependencies(baseline rrtmgp_allsky_baseline.nc)
