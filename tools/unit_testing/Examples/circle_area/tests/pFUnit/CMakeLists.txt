# Local pFUnit files.
set(pf_sources test_circle.pf)

# Sources to test.
set(sources_needed circle.F90)
extract_sources("${sources_needed}" "${circle_area_sources}" test_sources)

# Do source preprocessing and add the executable.
add_pFUnit_executable(pFUnit_circle_exe "${pf_sources}"
    ${CMAKE_CURRENT_BINARY_DIR} "${test_sources}")

add_test(pFUnit_circle_area pFUnit_circle_exe)

define_pFUnit_failure(pFUnit_circle_area)
