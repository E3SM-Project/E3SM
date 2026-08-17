set(SUPPORTS_CXX "TRUE")
#set(PIO_FILESYSTEM_HINTS "nfs") # need a check
#the following will override 'CMAKE_Fortran_FLAGS_DEBUG' from gnu.cmake, in which '-ffpt-trap=invalid' has seg. fault of 'mpi_init(ierr)' on pathfinder. A possible reason may be 'ierr' is 64-bit but should be 32-bit.
set(CMAKE_Fortran_FLAGS_DEBUG " -g -Wall -fbacktrace -fcheck=bounds,pointer -ffpe-trap=zero,overflow")
string(APPEND CMAKE_Fortran_FLAGS " -fallow-argument-mismatch")
#sometimes cannot find libMoab.so on pathfinder
set(MOAB_ROOT /projects/hpcl-cli185/proj-shared/ccsi-apps/moab)
if (COMP_INTERFACE STREQUAL "moab")
  if (COMP_NAME STREQUAL "cpl")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,-rpath,${MOAB_ROOT}/lib -L${MOAB_ROOT}/lib")
  endif()
endif()
