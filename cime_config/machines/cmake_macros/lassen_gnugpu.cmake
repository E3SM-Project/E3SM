string(APPEND CMAKE_C_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O2")
string(APPEND CMAKE_CUDA_FLAGS_RELEASE " -O3 -arch sm_70 --use_fast_math")
string(APPEND CMAKE_CUDA_FLAGS_DEBUG " -O0 -g -arch sm_70")

if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()

string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")

set(PIO_FILESYSTEM_HINTS "gpfs")

set(USE_CUDA "TRUE")

# This may not be needed once we figure out why MPI calls are segfaulting
# on lassen when this is ON.
set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "")
