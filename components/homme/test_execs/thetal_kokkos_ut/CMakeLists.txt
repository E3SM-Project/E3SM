SET(SRC_DIR          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR    ${HOMME_SOURCE_DIR}/src/share)
SET(SRC_THETA_DIR    ${HOMME_SOURCE_DIR}/src/theta-l_kokkos)
SET(THETA_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR})
SET(UTILS_SHARE_DIR  ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/cime/src/share/timing)

### ElementOps unit tests

SET (ELEM_OPS_UT_CXX_SRCS
  ${SRC_SHARE_DIR}/cxx/Context.cpp
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${THETA_UT_DIR}/elem_ops_ut.cpp
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
