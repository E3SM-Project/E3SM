set(test_sources test_infnan.F90)

set(sources_needed shr_kind_mod.F90 shr_infnan_mod.F90)

extract_sources("${sources_needed}" "${share_sources}" test_sources)

add_executable(infnan_test_exe ${test_sources})

declare_generated_dependencies(infnan_test_exe "${share_genf90_sources}")

if(NOT "${CMAKE_BUILD_TYPE}" MATCHES CESM_DEBUG)

  # Add the actual test.
  add_test(infnan infnan_test_exe)

  define_Fortran_stop_failure(infnan)

else()

  # CESM_DEBUG turns on floating point trapping, causing checks in this
  # test to fail. We'd like the users to be able to see that the test
  # was skipped but didn't fail, but CTest has no mechanism for this.
  # Instead, hack around it by adding a test that will always pass but
  # is likely to stand out as a skipped test.
  add_test(SKIPPED_infnan true)

endif()
