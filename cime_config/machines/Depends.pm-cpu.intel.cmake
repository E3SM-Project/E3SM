# For this file, we see internal compiler error with ifx (via intel-oneapi module) on pm-cpu with -O2
# Commenting for now as we are using intel module which is not seeing build issue
#set(NOOPT
#  eam/src/physics/cam/debug_info.F90)

#if (NOT DEBUG)
#  foreach(ITEM IN LISTS NOOPT)
#    e3sm_remove_flags("${ITEM}" "-O2")
#    e3sm_add_flags("${ITEM}" "-O0")
#  endforeach()
#endif()




