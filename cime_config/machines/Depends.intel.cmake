set(PERFOBJS
  homme/src/share/prim_advection_base.F90
  homme/src/share/vertremap_base.F90
  homme/src/share/edge_mod.F90
  homme/src/share/derivative_mod.F90
  homme/src/share/bndry_mod.F90
  homme/src/theta-l/share/prim_advance_mod.F90
  homme/src/preqx/share/prim_advance_mod.F90
  homme/src/preqx/share/viscosity_preqx_base.F90
  homme/src/share/viscosity_base.F90
  homme/src/theta-l/share/viscosity_theta.F90
  homme/src/theta-l/share/eos.F90
  eam/src/physics/cam/uwshcu.F90)

if (NOT DEBUG)
  foreach(ITEM IN LISTS PERFOBJS)
    e3sm_add_flags("${ITEM}" "-O3 -fp-model fast -no-prec-div")
  endforeach()
endif()

# https://github.com/E3SM-Project/E3SM/issues/8140
if (DEBUG)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
    e3sm_add_flags("eam/src/physics/cosp2/external/src/simulator/MISR_simulator/MISR_simulator.F90" "-init=zero")
  endif()
endif()

# https://github.com/E3SM-Project/E3SM/issues/8036
# Turn off vectorization in one source file to avoid issue with Intel compiler
if (NOT DEBUG)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "2025.0")
    e3sm_add_flags("eam/src/physics/cosp2/optics/quickbeam_optics.F90" "-no-vec")
  endif()
endif()
