# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  if (EAMXX_ENABLE_PYTHON)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/my_module.py
                    ${CMAKE_CURRENT_BINARY_DIR}/my_module.py)
    CreateUnitTest(pysession "pysession.cpp")
  endif()

  # Tests to ensure that fails/errors are indeed caught
  # NOTE: we don't need any libs for this test, but scream's CreateUnitTest
  # has libs as a required arg. So use the raw EkatCreateUnitTest

  # This serves as a base case for the following fail checks, since it verifies that
  #  - a REQUIRE clause that fails makes the test fail
  #  - Our Create unit test logic does work for  catching failures
  CreateUnitTest (fail "fail.cpp"
    WILL_FAIL LABELS "fail")

  if (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
    # Ensure that Kokkos OOB are caught
    CreateUnitTest (kokkos_fail "kokkos_fail.cpp"
      WILL_FAIL LABELS "fail")
  endif()

  if (EKAT_ENABLE_VALGRIND)
    # Ensure that valgrind errors are caught
    EkatCreateUnitTest (valg_fail "valg_fail.cpp"
      WILL_FAIL LABELS "fail")
  endif()

  # Ensure that FPE *do* throw when we expect them to
  if (SCREAM_FPE)
    CreateUnitTestExec (fpe_check "fpe_check.cpp"
      WILL_FAIL LABELS "check")
  else()
    CreateUnitTestExec (fpe_check "fpe_check.cpp"
      LABELS "check")
  endif()
endif()
