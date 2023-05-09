# For this file, see internal compiler error on pm-cpu with -O2
set(NOOPT
  eam/src/physics/cam/debug_info.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS NOOPT)
    e3sm_remove_flags("${ITEM}" "-O2")
    e3sm_add_flags("${ITEM}" "-O0")
  endforeach()
endif()

# compile mpas_seaice_core_interface.f90 with ifort, not ifx
e3sm_add_flags("${CMAKE_BINARY_DIR}/core_seaice/model_forward/mpas_seaice_core_interface.f90" "-fc=ifort")
