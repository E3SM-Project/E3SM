if (SCREAM_ONLY_GENERATE_BASELINES)
  # Generate allsky baseline
  # Note: these "baselines" are not to compare scream with a previous version, but
  #       rather to compare scream::rrtmgp with raw rrtmgp.
  CreateUnitTest (
    rrtmgp-allsky-baseline generate_baseline.cpp
    LIBS scream_rrtmgp rrtmgp_test_utils
    LABELS baseline_gen rrtmgp
    EXCLUDE_MAIN_CPP
    EXE_ARGS "${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc ${SCREAM_BASELINES_DIR}/data/rrtmgp-allsky-baseline.nc"
  )

else ()

  if (SCREAM_ENABLE_BASELINE_TESTS)
    # NOTE: tests inside this branch won't be built in a baselines-only build
    CreateUnitTest(rrtmgp_tests rrtmgp_tests.cpp
        LIBS scream_rrtmgp rrtmgp_test_utils
        LABELS "rrtmgp;physics"
        EXE_ARGS "-i ${SCREAM_DATA_DIR}/init/rrtmgp-allsky.nc -b ${SCREAM_BASELINES_DIR}/data/rrtmgp-allsky-baseline.nc"
        EXCLUDE_MAIN_CPP
    )
  endif()

  CreateUnitTest(rrtmgp_unit_tests rrtmgp_unit_tests.cpp
      LIBS scream_rrtmgp rrtmgp_test_utils
      LABELS "rrtmgp;physics"
  )
endif()
