# Reduce optimization on files with long compile times
set(O2MODELSRC
  eam/src/physics/cam/micro_mg_cam.F90      # ~2027 seconds
  elm/src/data_types/VegetationDataType.F90 # ~ 930 seconds
)
if (NOT DEBUG)
  foreach(ITEM IN LISTS O2MODELSRC)
    e3sm_remove_flags("${ITEM}" "-O3")
    e3sm_add_flags("${ITEM}" "-O2")
  endforeach()
endif()

