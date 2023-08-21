if (NOT DEBUG)
  string(APPEND CFLAGS " -O2")
  string(APPEND FFLAGS " -O2")
  string(APPEND CUDA_FLAGS " -O3 -arch sm_70 --use_fast_math")
endif()
if (DEBUG)
  string(APPEND CUDA_FLAGS " -O0 -g -arch sm_70")
endif()

if (COMP_NAME STREQUAL gptl)
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC")
endif()

string(APPEND CPPDEFS " -DTHRUST_IGNORE_CUB_VERSION_CHECK")
string(APPEND SLIBS " -L$ENV{NETCDF_PATH}/lib -lhdf5_hl -lhdf5 -lpnetcdf -lnetcdf -lnetcdff -lblas -llapack")

set(PIO_FILESYSTEM_HINTS "gpfs")

set(NETCDF_PATH "$ENV{NETCDF_PATH}")
set(NETCDF_C_PATH "$ENV{NETCDF_PATH}")
set(NETCDF_FORTRAN_PATH "$ENV{NETCDF_PATH}")
set(HDF5_PATH "$ENV{NETCDF_PATH}")
set(PNETCDF_PATH "$ENV{NETCDF_PATH}")
set(USE_CUDA "TRUE")

# This may not be needed once we figure out why MPI calls are segfaulting
# on lassen when this is ON.
set(SCREAM_MPI_ON_DEVICE OFF CACHE STRING "")
