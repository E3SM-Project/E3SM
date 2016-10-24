set(test_sources test_vmath.F90)

set(sources_needed shr_kind_mod.F90 shr_log_mod.F90 shr_strconvert_mod.F90 
      shr_infnan_mod.F90 shr_const_mod.F90 shr_vmath_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

add_executable(vmath_test_exe ${test_sources})

# Add the actual test.
add_test(vmath vmath_test_exe)

declare_generated_dependencies(vmath_test_exe "${share_genf90_sources}")
