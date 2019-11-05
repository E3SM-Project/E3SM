SET(SRC_DIR          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR    ${HOMME_SOURCE_DIR}/src/share)
SET(SRC_PREQX_DIR    ${HOMME_SOURCE_DIR}/src/preqx_kokkos)
SET(PREQX_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR})
SET(UTILS_SHARE_DIR  ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/cime/src/share/timing)

### Random initialization unit test

SET (RANDOM_INIT_UT_CXX_SRCS
  ${SRC_PREQX_DIR}/cxx/ElementsState.cpp
  ${SRC_SHARE_DIR}/cxx/Context.cpp
  ${SRC_SHARE_DIR}/cxx/ElementsGeometry.cpp
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${PREQX_UT_DIR}/random_init_ut.cpp
)
SET (RANDOM_INIT_UT_F90_SRCS)

SET (CONFIG_DEFINES PIO_INTERP NP=4 PLEV=8 QSIZE_D=4)
SET (RANDOM_INIT_UT_INCLUDE_DIRS
  ${SRC_PREQX_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${PREQX_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

SET (NUM_CPUS 1)
cxx_unit_test (random_init_ut "${RANDOM_INIT_UT_F90_SRCS}" "${RANDOM_INIT_UT_CXX_SRCS}" "${RANDOM_INIT_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})

### Remap unit test ###

SET (REMAP_PREQX_UT_F90_SRCS
  ${SRC_SHARE_DIR}/control_mod.F90
  ${SRC_SHARE_DIR}/dimensions_mod.F90
  ${SRC_SHARE_DIR}/hybvcoord_mod.F90
  ${SRC_SHARE_DIR}/kinds.F90
  ${SRC_SHARE_DIR}/parallel_mod.F90
  ${SRC_SHARE_DIR}/vertremap_base.F90
  ${SRC_SHARE_DIR}/physical_constants.F90
  ${UTILS_TIMING_DIR}/perf_mod.F90
  ${UTILS_TIMING_DIR}/perf_utils.F90
)

SET (REMAP_PREQX_UT_CXX_SRCS
  ${SRC_PREQX_DIR}/cxx/ElementsForcing.cpp
  ${SRC_PREQX_DIR}/cxx/ElementsState.cpp
  ${SRC_SHARE_DIR}/cxx/Context.cpp
  ${SRC_SHARE_DIR}/cxx/Elements.cpp
  ${SRC_SHARE_DIR}/cxx/ElementsDerivedState.cpp
  ${SRC_SHARE_DIR}/cxx/ElementsGeometry.cpp
  ${SRC_SHARE_DIR}/cxx/ErrorDefs.cpp
  ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  ${SRC_SHARE_DIR}/cxx/Hommexx_Session.cpp
  ${SRC_SHARE_DIR}/cxx/HybridVCoord.cpp
  ${SRC_SHARE_DIR}/cxx/Tracers.cpp
  ${SRC_SHARE_DIR}/cxx/mpi/Comm.cpp
  ${PREQX_UT_DIR}/remap_preqx_ut.cpp
)

SET (CONFIG_DEFINES NP=4 NC=4 PLEV=16 QSIZE_D=4 _MPI=1 _PRIM)
SET (REMAP_PREQX_UT_INCLUDE_DIRS
  ${SRC_PREQX_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${PREQX_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (remap_preqx_ut "${REMAP_PREQX_UT_F90_SRCS}" "${REMAP_PREQX_UT_CXX_SRCS}" "${REMAP_PREQX_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
