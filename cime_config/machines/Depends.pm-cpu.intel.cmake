# https://github.com/E3SM-Project/E3SM/issues/8140
if (DEBUG)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
    e3sm_add_flags("eam/src/physics/cosp2/external/src/simulator/MISR_simulator/MISR_simulator.F90" "-init=zero")
  endif()
endif()

# https://github.com/E3SM-Project/E3SM/issues/8036
if (NOT DEBUG)
  foreach(ITEM IN LISTS REDUCEDOPT)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
    e3sm_add_flags("eam/src/physics/cosp2/optics/quickbeam_optics.F90" "-O1")
  endif()
endif()
