# Reduce optimization on files with long compile times
# Also reduce optimizaiton on files with machine-dependent behavior
set(O2MODELSRC
  eam/src/physics/cam/micro_mg_cam.F90                            # ~ 2027 seconds
  elm/src/data_types/VegetationDataType.F90                       # ~ 930  seconds
  elm/src/external_models/fates/main/FatesHistoryInterfaceMod.F90 # ~ 1685 seconds
  elm/src/data_types/ColumnDataType.F90                           # ~ 556  seconds 
  eam/src/physics/cam/zm_conv.F90                                 # for machine dependent behavior
  eam/src/physics/cam/zm_microphysics.F90                         # for machine dependent behavior
)
if (NOT DEBUG)
  foreach(ITEM IN LISTS O2MODELSRC)
    e3sm_add_flags("${ITEM}" "-O2")
  endforeach()
endif()

