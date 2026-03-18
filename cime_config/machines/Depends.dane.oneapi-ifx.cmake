# compile mpas_seaice_core_interface.f90 with ifort, not ifx
if (NOT MPILIB STREQUAL "openmpi")
  e3sm_add_flags("${CMAKE_BINARY_DIR}/core_seaice/model_forward/mpas_seaice_core_interface.f90" "-fc=ifort")
  e3sm_add_flags("${CMAKE_BINARY_DIR}/core_landice/mode_forward/mpas_li_core_interface.f90" "-fc=ifort")
endif()
