include(EkatCreateUnitTest)

# Tests to ensure that fails/errors are indeed caught
# NOTE: we don't need any libs for this test, but scream's CreateUnitTest
# has libs as a required arg. So use the raw EkatCreateUnitTest

# Ensure that a non-satisfied REQUIRE clause does, in fact, make the test fail
EkatCreateUnitTestExec (fail "fail.cpp")
EkatCreateUnitTestFromExec (fail fail PROPERTIES WILL_FAIL TRUE)

if (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
  # Ensure that Kokkos OOB are caught
  EkatCreateUnitTestExec (kokkos_fail "kokkos_fail.cpp")
  EkatCreateUnitTestFromExec (kokkos_fail kokkos_fail PROPERTIES WILL_FAIL TRUE)
endif()

if (EKAT_ENABLE_VALGRIND)
  # Ensure that valgrind errors are caught
  EkatCreateUnitTestExec (valg_fail "valg_fail.cpp")
  EkatCreateUnitTestFromExec (valg_fail valg_fail PROPERTIES WILL_FAIL TRUE)
endif()

if (SCREAM_FPE)
  # Ensure that FPE *do* throw
  EkatCreateUnitTestExec (fpe_fail "fpe_fail.cpp")
  EkatCreateUnitTestFromExec (fpe_fail fpe_fail PROPERTIES WILL_FAIL TRUE)
endif()
