e3sm_add_flags("eam/src/dynamics/fv/geopk.F90" "-fcray-pointer")

list(APPEND MPAS_ICE_SHORTWAVE
  ${CMAKE_BINARY_DIR}/core_seaice/column/ice_shortwave.f90
)

if (NOT DEBUG)
  foreach(ITEM IN LISTS MPAS_ICE_SHORTWAVE)
    e3sm_add_flags("${ITEM}" "-fno-tree-pta") # avoids an error that shows up in solver with gnu9 and higher
  endforeach()
endif()
