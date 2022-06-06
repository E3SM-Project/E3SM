string(APPEND LDFLAGS " -gpu=cc70,cc60 -Minfo=accel")
string(APPEND SLIBS " -L$ENV{PNETCDF_PATH}/lib -lpnetcdf -L$ENV{HDF5_PATH}/lib -lhdf5_hl -lhdf5 -L$ENV{NETCDF_C_PATH}/lib -lnetcdf -L$ENV{NETCDF_FORTRAN_PATH}/lib -lnetcdff -L$ENV{ESSL_PATH}/lib64 -lessl -L$ENV{OLCF_NETLIB_LAPACK_ROOT}/lib64 -llapack")
set(KOKKOS_OPTIONS "--arch=Power9,Volta70 --with-cuda=$ENV{CUDA_DIR} --with-cuda-options=enable_lambda")
