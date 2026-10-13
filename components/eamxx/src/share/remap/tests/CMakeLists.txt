# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test iop remap
  CreateUnitTest(iop_remapper
    SOURCES iop_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test coarsening remap
  CreateUnitTest(coarsening_remapper
    SOURCES coarsening_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test refining remap
  CreateUnitTest(refining_remapper
    SOURCES refining_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test vertical remap
  CreateUnitTest(vertical_remapper
    SOURCES vertical_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
  )

  # Test remap to structured grids (lat-lon and curvilinear, e.g. HRRR-like)
  set(CURVILINEAR_MAP_FILE map_ne4pg2_to_10x20_hrrr_20261013.nc)
  set(CURVILINEAR_MAP_FILE_FULL ${SCREAM_INPUT_ROOT}/atm/scream/maps/${CURVILINEAR_MAP_FILE})
  GetInputFile ("scream/maps/${CURVILINEAR_MAP_FILE}")
  # Also test a lat-lon map file, since lat-lon and structured grids share the setup code
  set(LATLON_MAP_FILE map_ne4pg2_to_10x20_20260112.nc)
  set(LATLON_MAP_FILE_FULL ${SCREAM_INPUT_ROOT}/atm/scream/maps/${LATLON_MAP_FILE})
  GetInputFile ("scream/maps/${LATLON_MAP_FILE}")

  CreateUnitTest(horiz_structured_remap
    SOURCES horiz_structured_remap_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    EXE_ARGS --args --curvilinear-map-file=${CURVILINEAR_MAP_FILE_FULL} --latlon-map-file=${LATLON_MAP_FILE_FULL}
  )
endif()
