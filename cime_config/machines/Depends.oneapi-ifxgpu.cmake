
# compile mpas_seaice_core_interface.f90 with ifort, not ifx
if (NOT MPILIB STREQUAL "openmpi")
  e3sm_add_flags("${CMAKE_BINARY_DIR}/core_seaice/model_forward/mpas_seaice_core_interface.f90" "-fc=ifort")
endif()
