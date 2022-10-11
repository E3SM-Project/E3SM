# Some tests for checking testing works
set(EAMxx_tests_IC_FILE_72lev "screami_unit_tests_ne4np4L72_20220822.nc")
set(EAMxx_tests_IC_FILE_128lev "screami_unit_tests_ne2np4L128_20220822.nc")

add_subdirectory(generic)

if (NOT DEFINED ENV{SCREAM_FAKE_ONLY})
  # memcheck builds (and also coverage ones) can just run the max ranks, since they
  # do only need to perform checks on the code itself, rather than the model it represents
  if (EKAT_ENABLE_VALGRIND OR EKAT_ENABLE_CUDA_MEMCHECK OR EKAT_ENABLE_COVERAGE)
    set (TEST_RANK_START ${SCREAM_TEST_MAX_RANKS})
    set (TEST_RANK_END   ${SCREAM_TEST_MAX_RANKS})
  else()
    set (TEST_RANK_START 1)
    set (TEST_RANK_END   ${SCREAM_TEST_MAX_RANKS})
  endif()

  # Testing individual atm processes
  add_subdirectory(uncoupled)

  # Testing multiple atm processes coupled together
  add_subdirectory(coupled)
endif()
