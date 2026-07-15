# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci-snl.cmake)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-snl-oneapi CACHE STRING "")

# We use MKL instead of blas here
option(HOMME_USE_MKL "Whether to use Intel's MKL/oneMKL instead of blas/lapack" ON)

# Set kokkos arch
option (Kokkos_ARCH_SPR "" ON)

# Currently, we have 192 cores on blake's GraniteRapids nodes, but 4 ranks and 8 threads is plenty
set(SCREAM_TEST_MAX_RANKS   4 CACHE STRING "Upper limit on number of ranks for mpi tests")
set(SCREAM_TEST_MAX_THREADS 8 CACHE STRING "Upper limit on number of threads for threaded tests")
set(SCREAM_TEST_THREAD_INC  2 CACHE STRING "Thread count increment for threaded tests")

option (EAMXX_ENABLE_PYTHON "Whether to enable python interface from eamxx" ON)
