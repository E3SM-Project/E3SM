# For this file, see internal compiler error on pm-cpu with -O2
set(NOOPT
  eam/src/physics/cam/debug_info.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS NOOPT)
    e3sm_remove_flags("${ITEM}" "-O2")
    e3sm_add_flags("${ITEM}" "-O0")
  endforeach()
endif()




