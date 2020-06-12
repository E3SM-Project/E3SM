# Set defaults for EKAT testing
set (EKAT_TEST_MAX_THREADS 16 CACHE STRING "Max number of threads for tests")
set (EKAT_TEST_THREAD_INC 1  CACHE STRING "Increment for threads for tests")

# Process all subdirs
add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(kokkos)
add_subdirectory(pack)
add_subdirectory(units)
add_subdirectory(utils)
