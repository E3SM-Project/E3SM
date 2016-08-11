set(test_sources test_driver.F90)

extract_sources("interpolate_1d.F90" "${interpolate_sources}" test_sources)

add_executable(CTest_interpolate_exe ${test_sources})

add_test(interpolate CTest_interpolate_exe)

define_Fortran_stop_failure(interpolate)
