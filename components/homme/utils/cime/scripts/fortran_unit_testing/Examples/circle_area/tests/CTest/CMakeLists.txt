# Local files (just the test driver in this case).
set(test_sources test_driver.F90)

# Add sources we need from other directories (just circle.F90 in this
# case).
set(sources_needed circle.F90)
extract_sources("${sources_needed}" "${circle_area_sources}" test_sources)

add_executable(CTest_circle_exe ${test_sources})

add_test(circle_area CTest_circle_exe)

# Tell CTest how to figure out that "STOP 1" fails for the current
# compiler.
define_Fortran_stop_failure(circle_area)
