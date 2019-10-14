#Currently Loaded Modulefiles:
#  1) tcl/8.6.5
#  2) x11/libXau/1.0.7
#  3) x11/libX11/1.5.0
#  4) x11/xtrans/1.2.7
#  5) x11/xcb/1.8.1
#  6) x11/xproto/7.0.23
#  7) x11/xcb-proto/1.7.1
#  8) x11/pthread-stubs/0.3.0
#  9) tk/8.6.5
# 10) git/2.8.2
# 11) binutils/2.26.0
# 12) mpc/1.0.3
# 13) mpfr/3.1.3
# 14) gmp/6.1.0
# 15) isl/0.14.0
# 16) gdb/7.11.0
# 17) cmake/3.5.2
# 18) zlib/1.2.8
# 19) gcc/4.9.3
# 20) intel/compilers/17.2.174
# 21) java/oracle/1.7.0.75
# 22) numa/2.0.10
# 23) hwloc/1.11.4
# 24) openmpi/1.10.6/intel/17.2.174
# 25) hdf5/1.8.17/openmpi/1.10.4/intel/17.0.098
# 26) pnetcdf/1.7.0/openmpi/1.10.4/intel/17.0.098
# 27) netcdf/4.4.1/openmpi/1.10.4/intel/17.0.098

SET (CMAKE_BUILD_TYPE "" CACHE STRING "")
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

SET (ADD_Fortran_FLAGS "-xMIC_AVX512 -O3 -g" CACHE STRING "")
SET (ADD_C_FLAGS "-xMIC_AVX512 -O3 -g" CACHE STRING "")
SET (ADD_CXX_FLAGS "-xMIC_AVX512 -O3 -g" CACHE STRING "")

SET (ZLIB_DIR $ENV{ZLIB_ROOT} CACHE FILEPATH "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(USE_MPI_OPTIONS "--bind-to core" CACHE STRING "")


