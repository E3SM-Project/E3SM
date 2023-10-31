list(APPEND NOOPT_FILES eam/src/physics/cam/gw_common.F90)

set(O0MODELSRC
  eam/src/chemistry/aerosol/dust_sediment_mod.F90
  eam/src/chemistry/modal_aero/modal_aero_convproc.F90
  eam/src/chemistry/utils/modal_aero_calcsize.F90
  eam/src/physics/cam/zm_conv.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS O0MODELSRC)
    e3sm_remove_flags("${ITEM}" "-O1")
    e3sm_add_flags("${ITEM}" "-O0 -vector0")
  endforeach()

endif()

