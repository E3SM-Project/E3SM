SET(SRC_DIR          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR    ${HOMME_SOURCE_DIR}/src/share)
SET(SRC_THETA_DIR    ${HOMME_SOURCE_DIR}/src/theta-l_kokkos)
SET(SHARE_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR}/../share_kokkos_ut)
SET(THETA_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR})
SET(UTILS_SHARE_DIR  ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/cime/src/share/timing)

### Equation of state unit tests

SET (EOS_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/Context.cpp
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${THETA_UT_DIR}/eos_ut.cpp
)

SET (EOS_UT_F90_SRCS
  ${SRC_THETA_DIR}/element_state.F90
  ${SRC_THETA_DIR}/eos.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/edgetype_mod.F90
  ${SRC_SHARE_DIR}/element_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
  ${SRC_SHARE_DIR}/hybvcoord_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${THETA_UT_DIR}/eos_interface.F90
)

SET (CONFIG_DEFINES PIO_INTERP NP=4 PLEV=11 QSIZE_D=4 _MPI=1)
SET (EOS_UT_INCLUDE_DIRS
  ${SRC_THETA_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${THETA_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (eos_ut "${EOS_UT_F90_SRCS}" "${EOS_UT_CXX_SRCS}" "${EOS_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### Element ops unit tests

SET (ELEM_OPS_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/Context.cpp
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${THETA_UT_DIR}/elem_ops_ut.cpp
)

SET (ELEM_OPS_UT_F90_SRCS
  ${SRC_THETA_DIR}/element_state.F90
  ${SRC_THETA_DIR}/element_ops.F90
  ${SRC_THETA_DIR}/eos.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/edgetype_mod.F90
  ${SRC_SHARE_DIR}/element_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
  ${SRC_SHARE_DIR}/hybvcoord_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${UTILS_TIMING_DIR}/perf_mod.F90
  ${UTILS_TIMING_DIR}/perf_utils.F90
  ${THETA_UT_DIR}/elem_ops_interface.F90
)

SET (CONFIG_DEFINES PIO_INTERP NP=4 PLEV=11 QSIZE_D=4 _MPI=1)
SET (ELEM_OPS_UT_INCLUDE_DIRS
  ${SRC_THETA_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${THETA_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (elem_ops_ut "${ELEM_OPS_UT_F90_SRCS}" "${ELEM_OPS_UT_CXX_SRCS}" "${ELEM_OPS_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### HyperViscosity unit tests

SET (HYPERVISCOSITY_UT_CXX_SRCS
  ${SRC_THETA_DIR}/cxx/ElementsForcing.cpp
  ${SRC_THETA_DIR}/cxx/ElementsState.cpp
  ${SRC_THETA_DIR}/cxx/HyperviscosityFunctorImpl.cpp
  ${SRC_SHARE_DIR}/cxx/Context.cpp
  ${SRC_SHARE_DIR}/cxx/ElementsDerivedState.cpp
  ${SRC_SHARE_DIR}/cxx/ElementsGeometry.cpp
  ${SRC_SHARE_DIR}/cxx/Elements.cpp
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/FunctorsBuffersManager.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/HyperviscosityFunctor.cpp
  ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
  ${SRC_SHARE_DIR}/cxx/ReferenceElement.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/BoundaryExchange.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Connectivity.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/MpiBuffersManager.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/mpi_cxx_f90_interface.cpp
  ${THETA_UT_DIR}/hv_ut.cpp
)

SET (HYPERVISCOSITY_UT_F90_SRCS
  ${SRC_DIR}/repro_sum_mod.F90
  ${SRC_THETA_DIR}/bndry_mod.F90
  ${SRC_THETA_DIR}/derivative_mod.F90
  ${SRC_THETA_DIR}/edge_mod.F90
  ${SRC_THETA_DIR}/element_state.F90
  ${SRC_THETA_DIR}/viscosity_theta.F90
  ${SRC_SHARE_DIR}/bndry_mod_base.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/cube_mod.F90
  ${SRC_SHARE_DIR}/derivative_mod_base.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/edge_mod_base.F90
  ${SRC_SHARE_DIR}/edgetype_mod.F90
  ${SRC_SHARE_DIR}/element_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
  ${SRC_SHARE_DIR}/hybvcoord_mod.F90
  ${SRC_SHARE_DIR}/hybrid_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/ll_mod.F90
  ${SRC_SHARE_DIR}/mass_matrix_mod.F90
  ${SRC_SHARE_DIR}/metagraph_mod.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${SRC_SHARE_DIR}/quadrature_mod.F90
  ${SRC_SHARE_DIR}/schedtype_mod.F90
  ${SRC_SHARE_DIR}/schedule_mod.F90
  ${SRC_SHARE_DIR}/spacecurve_mod.F90
  ${SRC_SHARE_DIR}/thread_mod.F90
  ${SRC_SHARE_DIR}/utils_mod.F90
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90
  ${UTILS_TIMING_DIR}/perf_mod.F90
  ${UTILS_TIMING_DIR}/perf_utils.F90
  # ${SRC_THETA_DIR}/element_ops.F90
  # ${SRC_THETA_DIR}/eos.F90
  # ${THETA_UT_DIR}/elem_ops_interface.F90
  ${THETA_UT_DIR}/hv_interface.F90
  ${SHARE_UT_DIR}/geometry_interface.F90
)

SET (CONFIG_DEFINES PIO_INTERP NP=4 PLEV=11 QSIZE_D=4 _MPI=1)
SET (HYPERVISCOSITY_UT_INCLUDE_DIRS
  ${SRC_THETA_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${THETA_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (hv_ut "${HYPERVISCOSITY_UT_F90_SRCS}" "${HYPERVISCOSITY_UT_CXX_SRCS}" "${HYPERVISCOSITY_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
