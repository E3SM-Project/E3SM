# Files that cannot be compiled with -O2 without losing reproducibility
set(O1MODELSRC
  eam/src/chemistry/aerosol/dust_sediment_mod.F90
  eam/src/chemistry/modal_aero/modal_aero_convproc.F90
  eam/src/chemistry/utils/modal_aero_calcsize.F90
  eam/src/physics/cam/zm_conv.F90)

set(O0MODELSRC
  ../driver-mct/main/component_mod.F90
  ../driver-mct/main/cime_comp_mod.F90
)

if (NOT DEBUG)
  foreach(ITEM IN LISTS O1MODELSRC)
    e3sm_remove_flags("${ITEM}" "-O2")
    e3sm_add_flags("${ITEM}" "-O1 -Mnovect")
  endforeach()

  foreach(ITEM IN LISTS O0MODELSRC)
    e3sm_remove_flags("${ITEM}" "-O2")
    e3sm_remove_flags("${ITEM}" "-O1")
    e3sm_add_flags("${ITEM}" "-O0 -Mnovect")
  endforeach()
endif()
