SET(SRC_DIR          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR    ${HOMME_SOURCE_DIR}/src/share)
SET(SRC_PREQX_DIR    ${HOMME_SOURCE_DIR}/src/preqx_kokkos)
SET(PREQX_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR})
SET(UTILS_TIMING_DIR ${HOMME_SOURCE_DIR}/utils/cime/CIME/non_py/src/timing)

### Build a 'preqx_kokkos' library ###
# We do this once, and link against it for all unit tests. This reduces compilation times
PREQX_KOKKOS_SETUP()

# This is needed to compile the lib and test executables with the correct options
SET(THIS_CONFIG_IN ${HOMME_SOURCE_DIR}/src/preqx_kokkos/config.h.cmake.in)
SET(THIS_CONFIG_HC ${CMAKE_CURRENT_BINARY_DIR}/config.h.c)
SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)
SET (NUM_POINTS 4)
SET (NUM_PLEV 12)
SET (QSIZE_D 4)
SET (PIO_INTERP TRUE)
HommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

ADD_LIBRARY(preqx_kokkos_ut_lib
  ${PREQX_DEPS}
  ${TEST_SRC_F90}
  ${SRC_DIR}/checksum_mod.F90
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/common_movie_mod.F90
  ${SRC_DIR}/interpolate_driver_mod.F90
  ${SRC_DIR}/interp_movie_mod.F90
  ${SRC_DIR}/netcdf_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90
  ${SRC_DIR}/prim_movie_mod.F90
  ${SRC_DIR}/restart_io_mod.F90
  ${SRC_DIR}/surfaces_mod.F90
  ${SRC_DIR}/test_mod.F90
)
TARGET_INCLUDE_DIRECTORIES(preqx_kokkos_ut_lib PUBLIC ${EXEC_INCLUDE_DIRS})
TARGET_INCLUDE_DIRECTORIES(preqx_kokkos_ut_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_COMPILE_DEFINITIONS(preqx_kokkos_ut_lib PUBLIC "HAVE_CONFIG_H")
TARGET_LINK_LIBRARIES(preqx_kokkos_ut_lib Kokkos::kokkos)
TARGET_LINK_LIBRARIES(preqx_kokkos_ut_lib timing csm_share ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
IF (HOMME_USE_MKL)
  IF (MKL_TYPE STREQUAL "oneMKL")
    TARGET_LINK_LIBRARIES(preqx_kokkos_ut_lib -qmkl)
  ELSEIF (MKL_TYPE STREQUAL "Intel MKL")
    TARGET_LINK_LIBRARIES(preqx_kokkos_ut_lib -mkl)
  ENDIF ()
ENDIF()
IF(NOT BUILD_HOMME_WITHOUT_PIOLIBRARY)
  IF(HOMME_USE_SCORPIO)
    TARGET_LINK_LIBRARIES(preqx_kokkos_ut_lib piof pioc)
  ELSE ()
    TARGET_LINK_LIBRARIES(preqx_kokkos_ut_lib pio)
  ENDIF ()
ENDIF ()
# Fortran modules
SET(PREQX_LIB_MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/preqx_kokkos_ut_lib_modules)
SET_TARGET_PROPERTIES(preqx_kokkos_ut_lib PROPERTIES Fortran_MODULE_DIRECTORY ${PREQX_LIB_MODULE_DIR})

# Setting this once and for all, good for all unit tests.
SET (CONFIG_DEFINES HAVE_CONFIG_H)

### Random initialization unit test

SET (RANDOM_INIT_UT_CXX_SRCS
  ${PREQX_UT_DIR}/random_init_ut.cpp
)
SET (RANDOM_INIT_UT_F90_SRCS)

SET (RANDOM_INIT_UT_INCLUDE_DIRS
  ${SRC_PREQX_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${PREQX_UT_DIR}
  ${UTILS_TIMING_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (random_init_ut "${RANDOM_INIT_UT_F90_SRCS}" "${RANDOM_INIT_UT_CXX_SRCS}" "${RANDOM_INIT_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
TARGET_LINK_LIBRARIES(random_init_ut preqx_kokkos_ut_lib)

IF (HOMMEXX_BFB_TESTING)

  ### Remap unit test ###

  SET (REMAP_PREQX_UT_F90_SRCS
    # ${UTILS_TIMING_DIR}/perf_mod.F90
    # ${UTILS_TIMING_DIR}/perf_utils.F90
  )

  SET (REMAP_PREQX_UT_CXX_SRCS
    ${PREQX_UT_DIR}/remap_preqx_ut.cpp
  )

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
  TARGET_LINK_LIBRARIES(remap_preqx_ut preqx_kokkos_ut_lib)
ENDIF()
