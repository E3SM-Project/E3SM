list(APPEND NOOPT_FILES
  eam/src/physics/cam/micro_mg_cam.F90
  elm/src/data_types/ColumnDataType.F90
  elm/src/data_types/VegetationDataType.F90
  elm/src/biogeochem/CNNitrogenFluxType.F90
  elm/src/biogeochem/CNCarbonFluxType.F90
  mosart/src/wrm/WRM_subw_IO_mod.F90
  mosart/src/riverroute/RtmMod.F90
)

# Files added below to mitigate excessive compilation times
set(O1MODELSRC
  eam/src/physics/cam/micro_mg1_5.F90
  eam/src/physics/cam/micro_mg2_0.F90
)

set(O2MODELSRC
)

if (NOT DEBUG)
  foreach(ITEM IN LISTS O2MODELSRC)
    e3sm_add_flags("${ITEM}" "-O2")
  endforeach()
  foreach(ITEM IN LISTS O1MODELSRC)
    e3sm_add_flags("${ITEM}" "-O1")
  endforeach()
endif()


