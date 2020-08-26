# Files that cannot be compiled with -O2 without losing reproducibility
set(O1MODELSRC
  cam/src/chemistry/aerosol/dust_sediment_mod.F90
  cam/src/chemistry/modal_aero/modal_aero_convproc.F90
  cam/src/physics/cam/zm_conv.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS O1MODELSRC)
    e3sm_remove_flags("${ITEM}" "-O2")
    e3sm_add_flags("${ITEM}" "-O1 -Mnovect")
  endforeach()

endif()
