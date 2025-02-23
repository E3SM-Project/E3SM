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

list(APPEND NO_INLINE_ARG_PACKING_LIST
  eam/src/dynamics/se/inidat.F90
)

if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
  if (NOT DEBUG)
    # new in gnu10, inline arg packing was causing INF values with SMS_P4x1.ne4pg2_ne4pg2.F-MMFXX
    foreach(ITEM IN LISTS NO_INLINE_ARG_PACKING_LIST)
      e3sm_add_flags("${ITEM}" " -fno-inline-arg-packing")
    endforeach()
  endif()
endif()

