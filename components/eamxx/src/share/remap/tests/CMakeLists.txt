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

  # Test remap to a generic 2d rectilinear (HRRR-like) grid
  set(HRRR_MAP_FILE map_ne4pg2_to_hrrr_10x20_20261008.nc)
  set(HRRR_MAP_FILE_FULL ${SCREAM_INPUT_ROOT}/atm/scream/maps/${HRRR_MAP_FILE})
  GetInputFile ("scream/maps/${HRRR_MAP_FILE}")
  CreateUnitTest(rectilinear_remapper
    SOURCES rectilinear_remapper_tests.cpp
    LIBS eamxx_remap
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    EXE_ARGS --args --map-file=${HRRR_MAP_FILE_FULL} --nx=10 --ny=20
  )
endif()
