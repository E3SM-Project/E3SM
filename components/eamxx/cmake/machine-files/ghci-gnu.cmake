# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

# Set Fortran flags
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch" CACHE STRING "Fortran compiler flags" FORCE)

# Set the path to BLAS/LAPACK libraries
set(BLAS_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to BLAS library" FORCE)
set(LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/lib/libopenblas.so" CACHE STRING "Path to LAPACK library" FORCE)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-gnu CACHE STRING "")

set(EKAT_VALGRIND_SUPPRESSION_FILE "/projects/e3sm/data/baselines/scream/ghci-gnu/eamxx-valgrind.supp" CACHE FILEPATH "Use this valgrind suppression file if valgrind is enabled.")
option (EAMXX_ENABLE_PYTHON "Whether to enable python interface from eamxx" ON)
