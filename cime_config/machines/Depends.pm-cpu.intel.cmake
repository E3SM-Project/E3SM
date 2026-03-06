# For this file, we see internal compiler error with ifx (via intel-oneapi module) on pm-cpu with -O2
# Commenting for now as we are using intel module which is not seeing build issue
#set(NOOPT
#  eam/src/physics/cam/debug_info.F90)

#if (NOT DEBUG)
#  foreach(ITEM IN LISTS NOOPT)
#    e3sm_add_flags("${ITEM}" "-O0")
#  endforeach()
#endif()

# https://github.com/E3SM-Project/E3SM/issues/8036
set(INITZERO eam/src/physics/cosp2/external/src/simulator/MISR_simulator/MISR_simulator.F90)

if (DEBUG)
  foreach(ITEM IN LISTS INITZERO)
    e3sm_add_flags("${ITEM}" "-init=zero")
  endforeach()
endif()




