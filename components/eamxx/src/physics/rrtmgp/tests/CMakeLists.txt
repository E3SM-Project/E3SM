# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)

  # Required libraries
  set (NEED_LIBS scream_rrtmgp rrtmgp_test_utils)

  # Build baseline code
  add_executable(generate_baseline generate_baseline.cpp)
  target_link_libraries(generate_baseline PUBLIC ${NEED_LIBS})

  # Generate allsky baseline with the usual cmake custom command-target pair pattern
  # Note: these "baselines" are not to compare scream with a previous version, but
  #       rather to compare scream::rrtmgp with raw rrtmgp.
  add_custom_command (
      OUTPUT  ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc
      COMMAND ${CMAKE_COMMAND} -E env $<TARGET_FILE:generate_baseline>
              ${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc
              ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc
  )
  add_custom_target(rrtmgp_allsky_baseline.nc
      DEPENDS ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc
  )

  CreateUnitTest(
      rrtmgp_tests rrtmgp_tests.cpp "${NEED_LIBS}" LABELS "rrtmgp;physics"
      EXE_ARGS "-i ${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc -b ${SCREAM_TEST_DATA_DIR}/rrtmgp-allsky-baseline.nc"
      EXCLUDE_MAIN_CPP
  )
  add_dependencies (rrtmgp_tests rrtmgp_allsky_baseline.nc)

  CreateUnitTest(
      rrtmgp_unit_tests rrtmgp_unit_tests.cpp "${NEED_LIBS}" LABELS "rrtmgp;physics"
  )

endif()
