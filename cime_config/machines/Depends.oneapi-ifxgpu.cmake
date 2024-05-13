
set(CPPDEFS "${CPPDEFS} -DMPAS_OPENMP_OFFLOAD")
list(APPEND MPAS_ADD_ACC_FLAGS
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_mesh_pool.f90
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_velocity_solver_variational.f90
  ${CMAKE_BINARY_DIR}/core_seaice/shared/mpas_seaice_velocity_solver.f90
)

foreach(ITEM IN LISTS MPAS_ADD_ACC_FLAGS)
  e3sm_add_flags("${ITEM}" "-fiopenmp -fopenmp-targets=spir64")
endforeach()

# compile mpas_seaice_core_interface.f90 with ifort, not ifx
e3sm_add_flags("${CMAKE_BINARY_DIR}/core_seaice/model_forward/mpas_seaice_core_interface.f90" "-fc=ifort")
