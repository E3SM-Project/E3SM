if (compile_threaded)
  string(APPEND CMAKE_C_FLAGS   " -qopenmp")
  string(APPEND CMAKE_CXX_FLAGS " -qopenmp")
  string(APPEND CMAKE_Fortran_FLAGS   " -qopenmp")
  string(APPEND CMAKE_EXE_LINKER_FLAGS  " -qopenmp")
endif()

string(APPEND CMAKE_C_FLAGS " -fp-model precise -std=gnu99")
string(APPEND CMAKE_C_FLAGS_DEBUG   " -O0 -g")
string(APPEND CMAKE_C_FLAGS_RELEASE   " -O2 -debug minimal")

string(APPEND CMAKE_CXX_FLAGS_DEBUG " -O0 -g")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O2")

string(APPEND CMAKE_Fortran_FLAGS " -convert big_endian -assume byterecl -traceback -assume realloc_lhs")
string(APPEND CMAKE_Fortran_FLAGS_DEBUG   " -O0 -g ")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE   " -O2 -debug minimal")

string(APPEND CPPDEFS " -DFORTRANUNDERSCORE -DNO_R16 -DCPRINTEL")

string(APPEND CMAKE_Fortran_FORMAT_FIXED_FLAG " -fixed -132")
string(APPEND CMAKE_Fortran_FORMAT_FREE_FLAG " -free")
set(HAS_F2008_CONTIGUOUS "TRUE")
set(E3SM_LINK_WITH_FORTRAN "TRUE")

set(MPICC "mpicc")
set(MPICXX "mpicxx")
set(MPIFC "mpif90")

# Intel compilers have two different compiler families, the classic Intel 
# compilers (icc/icpc/ifort) and the newer Intel oneAPI compilers (icx/icpx/ifx).  
# The oneAPI compilers are based on LLVM and have a different compiler ID than 
# the classic Intel compilers.  Some flags for these two compiler families are 
# different, so we need to check which compiler family is being used and set the 
# flags accordingly.

# TODO: Remove these if statements once the classic Intel compilers are no longer supported.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  string(APPEND CMAKE_CXX_FLAGS " -fp-model precise")  
  string(APPEND CMAKE_Fortran_FLAGS " -fpscomp logicals -fp-model precise")
  string(APPEND CPPDEFS " -DHAVE_SLASHPROC -DHIDE_MPI")

  set(SCC "icx")
  set(SCXX "icpx")
  set(SFC "ifx")
else()
  string(APPEND CMAKE_CXX_FLAGS " -fp-model source")
  string(APPEND CMAKE_Fortran_FLAGS " -ftz -fp-model source")
  string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created")
  if (COMP_NAME STREQUAL cice)
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -init=nosnan,arrays")
  else()
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -init=snan,arrays")
  endif()

  set(SCC "icc")
  set(SCXX "icpc")
  set(SFC "ifort")
endif()
