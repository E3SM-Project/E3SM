set(SUPPORTS_CXX "TRUE")
#set(PIO_FILESYSTEM_HINTS "nfs") # need a check
string(APPEND CMAKE_Fortran_FLAGS " -fallow-argument-mismatch")
#sometimes cannot find libMoab.so on pathfinder
set(MOAB_ROOT /projects/hpcl-cli185/proj-shared/ccsi-apps/moab)
if (COMP_INTERFACE STREQUAL "moab")
  if (COMP_NAME STREQUAL "cpl")
    string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,-rpath,${MOAB_ROOT}/lib -L${MOAB_ROOT}/lib")
  endif()
endif()
