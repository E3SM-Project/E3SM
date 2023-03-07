include(ScreamUtils)

# Tests to ensure that fails/errors are indeed caught
# NOTE: we don't need any libs for this test, but scream's CreateUnitTest
# has libs as a required arg. So use the raw EkatCreateUnitTest

# Ensure that a non-satisfied REQUIRE clause does, in fact, make the test fail
EkatCreateUnitTestFromExec (fail "fail.cpp" PROPERTIES WILL_FAIL TRUE LABELS "fail")

if (Kokkos_ENABLE_DEBUG_BOUNDS_CHECK)
  # Ensure that Kokkos OOB are caught
  EkatCreateUnitTest (kokkos_fail "kokkos_fail.cpp" PROPERTIES WILL_FAIL TRUE LABELS "fail")
endif()

if (EKAT_ENABLE_VALGRIND)
  # Ensure that valgrind errors are caught
  EkatCreateUnitTest (valg_fail "valg_fail.cpp" PROPERTIES WILL_FAIL TRUE LABELS "fail")
endif()

# Ensure that FPE *do* throw when we expect them to
CreateUnitTestExec (fpe_check "fpe_check.cpp" "scream_share")
if (SCREAM_FPE)
  CreateUnitTestFromExec (fpe_check fpe_check PROPERTIES WILL_FAIL TRUE LABELS "check")
else()
  CreateUnitTestFromExec (fpe_check fpe_check LABELS "check")
endif()
