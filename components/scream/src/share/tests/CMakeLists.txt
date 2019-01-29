INCLUDE (ScreamUtils)

if (CUDA_BUILD)
  set(MAX_THREADS 1)
else()
  set(MAX_THREADS 16)
endif()

# Test utilities (f90)
CreateUnitTest(array_io "array_io_tests.cpp;unit_test.f90" scream_share)

# Test utilities (c++)
CreateUnitTest(util_cxx util_tests.cpp scream_share)

# Test packs
CreateUnitTest(packs "pack_tests.cpp;pack_kokkos_tests.cpp" scream_share)

# Test kokkos utils
CreateUnitTest(kokkos_utils "kokkos_utils_tests.cpp" scream_share THREADS 1 ${MAX_THREADS})

# Test fields
CreateUnitTest(field "field_tests.cpp" scream_share)

# Test workspace manager
CreateUnitTest(wsm "workspace_tests.cpp" scream_share THREADS 1 ${MAX_THREADS})
