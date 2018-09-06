SET(SRC_DIR          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR    ${HOMME_SOURCE_DIR}/src/share)
SET(SRC_PREQX_DIR    ${HOMME_SOURCE_DIR}/src/preqx_kokkos)
SET(SHARE_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR})
SET(UTILS_SHARE_DIR  ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/cime/src/share/timing)

### Infrastructure unit stest

SET (INFRASTRUCTURE_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${SHARE_UT_DIR}/infrastructure_ut.cpp
)
SET (INFRASTRUCTURE_UT_F90_SRCS)

SET (CONFIG_DEFINES PIO_INTERP NP=4 QSIZE_D=1 PLEV=8 HOMMEXX_NO_MPI_CONTEXT HOMMEXX_NO_CONTEXT)
SET (INFRASTRUCTURE_UT_INCLUDE_DIRS
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${SHARE_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

SET (NUM_CPUS 1)
cxx_unit_test (infrastructure_ut "${INFRASTRUCTURE_UT_F90_SRCS}" "${INFRASTRUCTURE_UT_CXX_SRCS}" "${INFRASTRUCTURE_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### Random initialization unit stest

SET (RANDOM_INIT_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Elements.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${SHARE_UT_DIR}/random_init_ut.cpp
)
SET (RANDOM_INIT_UT_F90_SRCS)

SET (CONFIG_DEFINES PIO_INTERP NP=4 PLEV=8 QSIZE_D=4 HOMMEXX_NO_CONTEXT HOMMEXX_NO_MPI_CONTEXT)
SET (RANDOM_INIT_UT_INCLUDE_DIRS
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${SHARE_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

SET (NUM_CPUS 1)
cxx_unit_test (random_init_ut "${RANDOM_INIT_UT_F90_SRCS}" "${RANDOM_INIT_UT_CXX_SRCS}" "${RANDOM_INIT_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### Boundary exchange unit test ###

SET (BOUNDARY_EXCHANGE_UT_F90_SRCS
  ${SRC_DIR}/repro_sum_mod.F90
  ${SRC_SHARE_DIR}/bndry_mod_base.F90
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/cube_mod.F90
  ${SRC_SHARE_DIR}/derivative_mod_base.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/edge_mod_base.F90
  ${SRC_SHARE_DIR}/edgetype_mod.F90
  ${SRC_SHARE_DIR}/element_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
  ${SRC_SHARE_DIR}/hybrid_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/ll_mod.F90
  ${SRC_SHARE_DIR}/metagraph_mod.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${SRC_SHARE_DIR}/quadrature_mod.F90
  ${SRC_SHARE_DIR}/schedtype_mod.F90
  ${SRC_SHARE_DIR}/schedule_mod.F90
  ${SRC_SHARE_DIR}/spacecurve_mod.F90
  ${SRC_SHARE_DIR}/thread_mod.F90
  ${SRC_SHARE_DIR}/viscosity_base.F90
  ${SRC_PREQX_DIR}/bndry_mod.F90
  ${SRC_PREQX_DIR}/derivative_mod.F90
  ${SRC_PREQX_DIR}/edge_mod.F90
  ${SRC_PREQX_DIR}/element_state.F90
  ${SHARE_UT_DIR}/boundary_exchange_ut.F90
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90
  ${UTILS_TIMING_DIR}/perf_mod.F90
  ${UTILS_TIMING_DIR}/perf_utils.F90
)
SET (BOUNDARY_EXCHANGE_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/mpi_cxx_f90_interface.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/BoundaryExchange.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/BuffersManager.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Connectivity.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/MpiContext.cpp
  ${SHARE_UT_DIR}/boundary_exchange_ut.cpp
)

SET (CONFIG_DEFINES NP=4 NC=4 PLEV=16 QSIZE_D=4 _MPI=1 _PRIM HOMMEXX_NO_CONTEXT)
SET (BOUNDARY_EXCHANGE_UT_INCLUDE_DIRS
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${SHARE_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (boundary_exchange_ut "${BOUNDARY_EXCHANGE_UT_F90_SRCS}" "${BOUNDARY_EXCHANGE_UT_CXX_SRCS}" "${BOUNDARY_EXCHANGE_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### Sphere operators unit test ###

SET (SPHERE_OP_UT_F90_SRCS
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/derivative_mod_base.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/edgetype_mod.F90
  ${SRC_SHARE_DIR}/element_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/quadrature_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${SRC_PREQX_DIR}/element_state.F90
  ${SHARE_UT_DIR}/sphere_op_interface.F90
)
SET (SPHERE_OP_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${SHARE_UT_DIR}/sphere_op_sl.cpp
  ${SHARE_UT_DIR}/sphere_op_ml.cpp
)

SET (CONFIG_DEFINES NP=4 NC=4 PLEV=16 QSIZE_D=4 _MPI=1 _PRIM HOMMEXX_NO_CONTEXT HOMMEXX_NO_MPI_CONTEXT)
SET (SPHERE_OP_UT_INCLUDE_DIRS
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${SHARE_UT_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (sphere_op_ut "${SPHERE_OP_UT_F90_SRCS}" "${SPHERE_OP_UT_CXX_SRCS}" "${SPHERE_OP_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### Limiters unit test ###

SET (LIMITERS_UT_F90_SRCS
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/coordinate_systems_mod.F90
  ${SRC_SHARE_DIR}/derivative_mod_base.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/edgetype_mod.F90
  ${SRC_SHARE_DIR}/element_mod.F90
  ${SRC_SHARE_DIR}/gridgraph_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${SRC_SHARE_DIR}/params_mod.F90
  ${SRC_SHARE_DIR}/quadrature_mod.F90
  ${SRC_PREQX_DIR}/element_state.F90
  ${SHARE_UT_DIR}/limiters_interface.F90
)
SET (LIMITERS_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/EulerStepFunctorImpl.hpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SHARE_UT_DIR}/limiters.cpp
)

SET (CONFIG_DEFINES NP=4 NC=4 PLEV=13 QSIZE_D=4 _MPI=1 _PRIM HOMMEXX_NO_CONTEXT HOMMEXX_NO_MPI_CONTEXT)
SET (LIMITERS_UT_INCLUDE_DIRS
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${SHARE_UT_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (limiters_ut "${LIMITERS_UT_F90_SRCS}" "${LIMITERS_UT_CXX_SRCS}" "${LIMITERS_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### Remap unit test ###

SET (REMAP_UT_F90_SRCS
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/hybvcoord_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/vertremap_base.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${SHARE_UT_DIR}/remap_interface.F90
  ${UTILS_TIMING_DIR}/perf_mod.F90
  ${UTILS_TIMING_DIR}/perf_utils.F90
)

SET (REMAP_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/Elements.cpp
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
  ${SRC_SHARE_DIR}/cxx/Tracers.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${SHARE_UT_DIR}/remap_ut.cpp
)

SET (CONFIG_DEFINES NP=4 NC=4 PLEV=16 QSIZE_D=4 _MPI=1 _PRIM HOMMEXX_NO_CONTEXT HOMMEXX_NO_MPI_CONTEXT)
SET (REMAP_UT_INCLUDE_DIRS
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${SHARE_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (remap_ut "${REMAP_UT_F90_SRCS}" "${REMAP_UT_CXX_SRCS}" "${REMAP_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
