# https://github.com/E3SM-Project/E3SM/issues/8036
set(INITZERO eam/src/physics/cosp2/external/src/simulator/MISR_simulator/MISR_simulator.F90)

if (DEBUG)
  foreach(ITEM IN LISTS INITZERO)
    e3sm_add_flags("${ITEM}" "-init=zero")
  endforeach()
endif()




