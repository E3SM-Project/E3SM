list(APPEND MPAS_ICE_SHORTWAVE
  ${CMAKE_BINARY_DIR}/core_seaice/column/ice_shortwave.f90
)

# For optimized GNU builds that use v9 or higher, remove an optimization on one file
if (NOT DEBUG)
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 9)
    foreach(ITEM IN LISTS MPAS_ICE_SHORTWAVE)
      e3sm_add_flags("${ITEM}" "-fno-tree-pta") # avoids an error that shows up in solver with gnu9 and higher
    endforeach()
  endif()
endif()

list(APPEND ALLOW_INVALID_BOZ_LIST
  eam/src/control/cam_history.F90
  elm/src/biogeochem/MEGANFactorsMod.F90
)

if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  foreach(ITEM IN LISTS ALLOW_INVALID_BOZ_LIST)
    e3sm_add_flags("${ITEM}" "-fallow-invalid-boz") # avoids build error for integer, parameter :: gen_hash_key_offset = z'000053db'
  endforeach()
endif()

