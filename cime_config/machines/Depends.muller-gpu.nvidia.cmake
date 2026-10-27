list(APPEND REDUCE_OPT_LIST
  homme/src/share/derivative_mod.F90
)

# -Mnovect avoids an internal compiler error in this file (with nvidia/21.11).
# NVHPC >= 25.9 has fixed the ICE; skip the flag for those versions.
if (NOT DEBUG AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "25.9")
  foreach(ITEM IN LISTS REDUCE_OPT_LIST)
    e3sm_add_flags("${ITEM}" " -Mnovect")
  endforeach()
endif()
