string(APPEND FFLAGS " -fp-model consistent -fimf-use-svml")
if (NOT DEBUG)
  string(APPEND FFLAGS " -qno-opt-dynamic-align")
endif()
string(APPEND CXXFLAGS " -fp-model consistent")
set(PETSC_PATH "$ENV{PETSC_DIR}")
set(SCC "icc")
set(SCXX "icpc")
set(SFC "ifort")

set(HDF5_PATH, "$ENV{HDF5_PATH}")
set(NETCDF_C_PATH, "$ENV{NETCDF_C_PATH}")
set(NETCDF_FORTRAN_PATH, "$ENV{NETCDF_FORTRAN_PATH}")
set(PNETCDF_PATH, "$ENV{PNETCDF_PATH}")
#string(APPEND SLIBS " -L$ENV{NETCDF_DIR} -lnetcdff -Wl,--as-needed,-L$ENV{NETCDF_DIR}/lib -lnetcdff -lnetcdf")
string(APPEND SLIBS " -L$ENV{HDF5_PATH}/lib -lhdf5_fortran -lhdf5 -lhdf5_hl -lhdf5hl_fortran")
string(APPEND SLIBS " -L$ENV{NETCDF_C_PATH}/lib -lnetcdf -L$ENV{NETCDF_FORTRAN_PATH}/lib -lnetcdff")
string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf")
string(APPEND SLIBS " -mkl")

#   <HDF5_PATH>$ENV{HDF5_PATH}</HDF5_PATH>
#   <NETCDF_C_PATH>$ENV{NETCDF_C_PATH}</NETCDF_C_PATH>
#   <NETCDF_FORTRAN_PATH>$ENV{NETCDF_FORTRAN_PATH}</NETCDF_FORTRAN_PATH>
#   <PNETCDF_PATH>$ENV{PNETCDF_PATH}</PNETCDF_PATH>
#   <MPICC MPILIB="impi"> mpiicc  </MPICC>
#   <MPICXX MPILIB="impi"> mpiicpc </MPICXX>
#   <MPIFC MPILIB="impi"> mpiifort </MPIFC>
