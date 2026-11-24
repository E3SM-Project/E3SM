# Unit tests of the C++ code. The utilities to create them are in UnitTestUtils.cmake.

add_subdirectory (share)

if (BUILD_HOMME_THETA_KOKKOS)
  add_subdirectory (theta-l_kokkos)
endif()

if (BUILD_HOMME_PREQX_KOKKOS AND HOMMEXX_BFB_TESTING)
  add_subdirectory (preqx_kokkos)
endif()
