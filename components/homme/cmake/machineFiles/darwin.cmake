# CMake initial cache file for Mac OSX 10.8
# tested with gcc/gfortran & openmpi from HOMEBREW
#
# openmpi by default will use Mac's clang based compiler. to use GCC:
#
# setenv OMPI_CXX g++-9
# setenv OMPI_CC gcc-9
# setenv OMPI_FC gfortran-9
#
# also of interest:
# setenv  OMPI_MCA_btl self,tcp    # workaround https://github.com/open-mpi/ompi/issues/6518
#

# openmp is supported by GNU gcc, but not supported by Apple LLVM gcc
SET (ENABLE_OPENMP FALSE CACHE FILEPATH "")

SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")

SET (ADD_Fortran_FLAGS "-fbacktrace" CACHE STRING "")
#SET (ADD_Fortran_FLAGS "-fbacktrace -fbounds-check -O0" CACHE STRING "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (NETCDF_DIR /usr/local CACHE FILEPATH "")
SET (HDF5_DIR /usr/local/ CACHE FILEPATH "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

# see ../../README.cmake for options regarding cprnc
# For non-CIME supported systems, assume cprnc is already built:
#SET (CPRNC_DIR /usr/local/bin CACHE FILEPATH "")

SET (USE_MPI_OPTIONS "-oversubscribe" CACHE FILEPATH "")