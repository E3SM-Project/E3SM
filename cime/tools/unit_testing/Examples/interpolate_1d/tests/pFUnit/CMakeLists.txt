set(pf_sources test_interpolate_point.pf)

extract_sources("interpolate_1d.F90" "${interpolate_sources}" test_sources)

add_pFUnit_executable(pFUnit_interpolate_exe "${pf_sources}"
    ${CMAKE_CURRENT_BINARY_DIR} "${test_sources}")

add_test(pFUnit_interpolate pFUnit_interpolate_exe)

define_pFUnit_failure(pFUnit_interpolate)
