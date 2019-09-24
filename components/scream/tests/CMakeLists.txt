IF ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
  add_subdirectory(scream_homme_dyn_ut_nlev72_qsize4)
  add_subdirectory(scream_homme_dyn_p3)
ENDIF()
# TODO: Fix the P3 tests that go through the scream-ad infrastructure.
#       They break when kokkos is set to a debug build.
#add_subdirectory(scream_p3)
#add_subdirectory(scream_p3_shoc)  
