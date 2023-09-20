# For this file, fixes non-BFB behavior of stealth feature on pm-cpu with -O2
set(NOOPT
  eam/src/physics/cam/zm_conv.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS NOOPT)
    e3sm_remove_flags("${ITEM}" "-O2")
    e3sm_remove_flags("${ITEM}" "-O")
  endforeach()
endif()




