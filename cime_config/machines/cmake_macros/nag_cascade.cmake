string(APPEND CPPDEFS " -DnoI8")
if (DEBUG)
  string(APPEND FFLAGS " -C=all  -g  -O0 -v")
endif()
if (DEBUG AND COMP_NAME STREQUAL eam)
  string(APPEND FFLAGS " -C=all  -g  -nan -O0 -v")
endif()
if (MPILIB STREQUAL mvapich2)
  set(MPI_PATH "$ENV{MPI_LIB}")
endif()
set(NETCDF_PATH "$ENV{NETCDF_ROOT}")
set(PIO_FILESYSTEM_HINTS "lustre")
set(PNETCDF_PATH "$ENV{PNETCDFROOT}")
string(APPEND SLIBS " -L$ENV{NETCDF_ROOT}/lib -lnetcdf -lnetcdff -L$ENV{MKL_PATH} -lmkl_rt")