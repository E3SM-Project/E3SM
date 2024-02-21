include (ScreamUtils)

# Some tests for checking that certain testing infrastructures work
add_subdirectory(generic)

if (NOT DEFINED ENV{SCREAM_FAKE_ONLY})
  # memcheck builds (and also coverage ones) can just run the max ranks, since they
  # do only need to perform checks on the code itself, rather than the model it represents
  SetVarDependingOnTestSize(TEST_RANK_START ${SCREAM_TEST_MAX_RANKS} 1 1)
  set(TEST_RANK_END ${SCREAM_TEST_MAX_RANKS})

  # Initial condition files used in the tests
  set(EAMxx_tests_IC_FILE_72lev  "screami_unit_tests_ne2np4L72_20220822.nc")
  set(EAMxx_tests_IC_FILE_128lev "screami_unit_tests_ne2np4L128_20220822.nc")
  set(EAMxx_tests_IC_FILE_MAM    "scream_unit_tests_aerosol_optics_ne2np4L72_20220822.nc")
  set(EAMxx_tests_TOPO_FILE      "USGS-gtopo30_ne2np4pg2_x6t_20230331.nc")

  # Testing individual atm processes
  add_subdirectory(uncoupled)

  # Testing multiple atm processes coupled together
  add_subdirectory(coupled)
endif()
