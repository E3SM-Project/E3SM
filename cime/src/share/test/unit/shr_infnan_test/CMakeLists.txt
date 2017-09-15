set(test_sources test_infnan.F90)

set(sources_needed shr_kind_mod.F90 shr_infnan_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

add_executable(infnan_test_exe ${test_sources})

declare_generated_dependencies(infnan_test_exe "${share_genf90_sources}")

# Add the actual test.
add_test(infnan infnan_test_exe)

define_Fortran_stop_failure(infnan)
