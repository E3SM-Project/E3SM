# ANL Mira/Cetus machines

SET (CMAKE_Fortran_COMPILER mpixlf90_r CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpixlc_r CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpixlcxx_r CACHE FILEPATH "")

SET (CMAKE_Fortran_FLAGS "-WF,-C! -O2 -Ipreqx_modules" CACHE STRING "")
#SET (FORCE_Fortran_FLAGS "-WF,-C!" CACHE STRING "")
SET (ENABLE_OPENMP TRUE CACHE BOOL "")

#SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (PNETCDF_DIR /soft/libraries/pnetcdf/current/cnk-xl/current CACHE FILEPATH "")
SET (NETCDF_DIR /soft/libraries/netcdf/current/cnk-xl/current CACHE FILEPATH "")
SET (HDF5_DIR /soft/libraries/hdf5/current/cnk-xl/current CACHE FILEPATH "")
SET (ZLIB_DIR /soft/libraries/alcf/current/xl/ZLIB CACHE FILEPATH "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (USE_MPIEXEC runjob CACHE FILEPATH "")

SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
SET (ENV{PATH} "/soft/libraries/alcf/current/xl/BLAS/lib:/soft/libraries/alcf/current/xl/LAPACK/lib:$ENV{PATH}")
SET (ENV{LD_LIBRARY_PATH} "/soft/libraries/alcf/current/xl/BLAS/lib:/soft/libraries/alcf/current/xl/LAPACK/lib:$ENV{LD_LIBRARY_PATH}")

