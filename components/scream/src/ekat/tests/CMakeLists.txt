include(ScreamUtils)

# Test debug tools
CreateUnitTest(debug_tools "debug_tools_tests.cpp" ekat)

# Test utilities (f90)
CreateUnitTest(array_io "array_io_tests.cpp;unit_test.f90" ekat)

# Test utilities (c++)
CreateUnitTest(util_cxx "util_tests.cpp" ekat)

# Test packs
CreateUnitTest(packs "pack_tests.cpp;pack_kokkos_tests.cpp" ekat)

# Test kokkos utils
CreateUnitTest(kokkos_utils "kokkos_utils_tests.cpp" ekat THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})

# Test workspace manager
CreateUnitTest(wsm "workspace_tests.cpp" ekat THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})

# Test lin interp processes
CreateUnitTest(lin_interp "lin_interp_test.cpp;lin_interp_ref.f90" ekat THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})

# Test upper_bound
CreateUnitTest(upper_bound "upper_bound_test.cpp" ekat)

# Test tridiag solvers
CreateUnitTest(tridiag "tridiag_tests.cpp;tridiag_tests_correctness.cpp;tridiag_tests_performance.cpp" ekat THREADS ${SCREAM_TEST_MAX_THREADS} OPTIONAL EXCLUDE_MAIN_CPP)

# Test units framework
CreateUnitTest(units "units.cpp" ekat)

# YAML parser
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)
CreateUnitTest(yaml_parser "yaml_parser.cpp" ekat)
