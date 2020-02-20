# anlworkstation

SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

#SET (CMAKE_Fortran_FLAGS "-WF,-C! -O2 -Ipreqx_modules" CACHE STRING "")
#SET (FORCE_Fortran_FLAGS "-WF,-C!" CACHE STRING "")
SET (ENABLE_OPENMP TRUE CACHE BOOL "")

SET (PNETCDF_DIR /soft/apps/packages/climate/pnetcdf/1.12.0/mpich-3.3.2/gcc-8.2.0 CACHE FILEPATH "")
SET (NETCDF_DIR /soft/apps/packages/climate/netcdf/4.4.1c-4.2cxx-4.4.4f-parallel/mpich-3.3.2/gcc-8.2.0 CACHE FILEPATH "")
SET (ENV{HDF5} /soft/apps/packages/climate/hdf5/1.8.16-parallel/mpich-3.3.2/gcc-8.2.0 CACHE FILEPATH "")
#SET (ENV{LIBZ} /soft/libraries/alcf/current/xl/ZLIB CACHE FILEPATH "")
SET (CPRNC_DIR /home/climate1/acme/cprnc/build/cprnc CACHE FILEPATH "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (USE_MPIEXEC mpiexec CACHE FILEPATH "")

SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
#SET (ENV{PATH} "/soft/libraries/alcf/current/xl/BLAS/lib:/soft/libraries/alcf/current/xl/LAPACK/lib:$ENV{PATH}")
SET (ENV{PATH} "/soft/apps/packages/climate/netcdf/4.4.1c-4.2cxx-4.4.4f-parallel/mpich-3.3.2:$ENV{PATH}")
SET (ENV{LD_LIBRARY_PATH} "/soft/apps/packages/climate/hdf5/1.8.16-serial/gcc-8.2.0/lib:$ENV{LD_LIBRARY_PATH}")

SET (BUILD_HOMME_SWEQX FALSE CACHE BOOL "")
