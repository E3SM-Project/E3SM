SET(SRC_DIR          ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR    ${HOMME_SOURCE_DIR}/src/share)
SET(SRC_THETA_DIR    ${HOMME_SOURCE_DIR}/src/theta-l_kokkos)
SET(SHARE_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR}/../share_kokkos_ut)
SET(THETA_UT_DIR     ${CMAKE_CURRENT_SOURCE_DIR})
SET(UTILS_TIMING_SRC_DIR ${HOMME_SOURCE_DIR}/utils/cime/CIME/non_py/src/timing)
SET(UTILS_TIMING_BIN_DIR ${HOMME_BINARY_DIR}/utils/cime/CIME/non_py/src/timing)

### Build a 'theta-l_kokkos' library ###
# We do this once, and link against it for all unit tests. This reduces compilation times
THETAL_KOKKOS_SETUP()

# This is needed to compile the lib and test executables with the correct options
#these vars shared between all targets, so changing one var
#for one test only won't work, config is built once and for the last test
SET(THIS_CONFIG_IN ${HOMME_SOURCE_DIR}/src/theta-l_kokkos/config.h.cmake.in)
SET(THIS_CONFIG_HC ${CMAKE_CURRENT_BINARY_DIR}/config.h.c)
SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)
SET (NUM_POINTS 4)
SET (NUM_PLEV 12)
SET (QSIZE_D 4)
SET (PIO_INTERP TRUE)

HommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

ADD_LIBRARY(thetal_kokkos_ut_lib
  ${THETAL_DEPS}
  ${TEST_SRC_F90}
  ${SRC_DIR}/checksum_mod.F90
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/common_movie_mod.F90
  ${SRC_DIR}/interpolate_driver_mod.F90
  ${SRC_DIR}/interp_movie_mod.F90
  ${SRC_DIR}/netcdf_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90
  ${SRC_DIR}/prim_movie_mod.F90
  ${SRC_DIR}/theta_restart_mod.F90
  ${SRC_DIR}/restart_io_mod.F90
  ${SRC_DIR}/surfaces_mod.F90
  ${SRC_DIR}/test_mod.F90
)
TARGET_INCLUDE_DIRECTORIES(thetal_kokkos_ut_lib PUBLIC ${EXEC_INCLUDE_DIRS})
TARGET_INCLUDE_DIRECTORIES(thetal_kokkos_ut_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
TARGET_COMPILE_DEFINITIONS(thetal_kokkos_ut_lib PUBLIC "HAVE_CONFIG_H")
TARGET_LINK_LIBRARIES(thetal_kokkos_ut_lib Kokkos::kokkos)
TARGET_LINK_LIBRARIES(thetal_kokkos_ut_lib timing csm_share ${COMPOSE_LIBRARY_CPP} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
IF (HOMME_USE_MKL)
  IF (MKL_TYPE STREQUAL "oneMKL")
    TARGET_LINK_LIBRARIES(thetal_kokkos_ut_lib -qmkl)
  ELSEIF (MKL_TYPE STREQUAL "Intel MKL")
    TARGET_LINK_LIBRARIES(thetal_kokkos_ut_lib -mkl)
  ENDIF ()
ENDIF()
IF(BUILD_HOMME_WITHOUT_PIOLIBRARY)
  TARGET_COMPILE_DEFINITIONS(thetal_kokkos_ut_lib PUBLIC HOMME_WITHOUT_PIOLIBRARY)
ELSE ()
  IF(HOMME_USE_SCORPIO)
    TARGET_LINK_LIBRARIES(thetal_kokkos_ut_lib piof pioc)
  ELSE ()
    TARGET_LINK_LIBRARIES(thetal_kokkos_ut_lib pio)
  ENDIF ()
ENDIF ()
# Fortran modules
SET(THETA_LIB_MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/thetal_kokkos_ut_lib_modules)
SET_TARGET_PROPERTIES(thetal_kokkos_ut_lib PROPERTIES Fortran_MODULE_DIRECTORY ${THETA_LIB_MODULE_DIR})

# Setting this once and for all, good for all unit tests.
SET (CONFIG_DEFINES HAVE_CONFIG_H)

IF (HOMMEXX_BFB_TESTING)
  ### Equation of state unit tests

  SET (EOS_UT_CXX_SRCS
    ${THETA_UT_DIR}/eos_ut.cpp
    )

  SET (EOS_UT_F90_SRCS
    ${THETA_UT_DIR}/eos_interface.F90
    )

  SET (EOS_UT_INCLUDE_DIRS
    ${SRC_THETA_DIR}/cxx
    ${SRC_SHARE_DIR}
    ${SRC_SHARE_DIR}/cxx
    ${THETA_UT_DIR}
    ${THETA_LIB_MODULE_DIR}
    ${UTILS_TIMING_SRC_DIR}
    ${UTILS_TIMING_BIN_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/src/share/cxx
    )

  IF (USE_NUM_PROCS)
    SET (NUM_CPUS ${USE_NUM_PROCS})
  ELSE()
    SET (NUM_CPUS 1)
  ENDIF()
  cxx_unit_test (eos_ut "${EOS_UT_F90_SRCS}" "${EOS_UT_CXX_SRCS}" "${EOS_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
  TARGET_LINK_LIBRARIES(eos_ut thetal_kokkos_ut_lib)

  ### Element ops unit tests

  SET (ELEM_OPS_UT_CXX_SRCS
    ${THETA_UT_DIR}/elem_ops_ut.cpp
    )

  SET (ELEM_OPS_UT_F90_SRCS
    ${THETA_UT_DIR}/elem_ops_interface.F90
    )

  SET (ELEM_OPS_UT_INCLUDE_DIRS
    ${SRC_THETA_DIR}/cxx
    ${SRC_SHARE_DIR}
    ${SRC_SHARE_DIR}/cxx
    ${THETA_UT_DIR}
    ${THETA_LIB_MODULE_DIR}
    ${UTILS_TIMING_SRC_DIR}
    ${UTILS_TIMING_BIN_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/src/share/cxx
    )

  IF (USE_NUM_PROCS)
    SET (NUM_CPUS ${USE_NUM_PROCS})
  ELSE()
    SET (NUM_CPUS 1)
  ENDIF()
  cxx_unit_test (elem_ops_ut "${ELEM_OPS_UT_F90_SRCS}" "${ELEM_OPS_UT_CXX_SRCS}" "${ELEM_OPS_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
  TARGET_LINK_LIBRARIES(elem_ops_ut thetal_kokkos_ut_lib)

  # ### HyperViscosity unit tests

  SET (HYPERVISCOSITY_UT_CXX_SRCS
    ${THETA_UT_DIR}/hv_ut.cpp
    )

  SET (HYPERVISCOSITY_UT_F90_SRCS
    ${THETA_UT_DIR}/hv_interface.F90
    ${THETA_UT_DIR}/thetal_test_interface.F90
    ${SHARE_UT_DIR}/geometry_interface.F90
    )

  SET (HYPERVISCOSITY_UT_INCLUDE_DIRS
    ${SRC_THETA_DIR}/cxx
    ${SRC_SHARE_DIR}
    ${SRC_SHARE_DIR}/cxx
    ${THETA_UT_DIR}
    ${THETA_LIB_MODULE_DIR}
    ${UTILS_TIMING_SRC_DIR}
    ${UTILS_TIMING_BIN_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/src/share/cxx
    )

  IF (USE_NUM_PROCS)
    SET (NUM_CPUS ${USE_NUM_PROCS})
  ELSE()
    SET (NUM_CPUS 1)
  ENDIF()
  cxx_unit_test (hv_ut "${HYPERVISCOSITY_UT_F90_SRCS}" "${HYPERVISCOSITY_UT_CXX_SRCS}" "${HYPERVISCOSITY_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
  TARGET_LINK_LIBRARIES(hv_ut thetal_kokkos_ut_lib)

  ### Forcing unit tests

  SET (FORCING_UT_CXX_SRCS
    ${THETA_UT_DIR}/forcing_ut.cpp
    )

  SET (FORCING_UT_F90_SRCS
    ${THETA_UT_DIR}/forcing_interface.F90
    ${SHARE_UT_DIR}/geometry_interface.F90
    )

  SET (FORCING_UT_INCLUDE_DIRS
    ${SRC_THETA_DIR}/cxx
    ${SRC_SHARE_DIR}
    ${SRC_SHARE_DIR}/cxx
    ${THETA_UT_DIR}
    ${THETA_LIB_MODULE_DIR}
    ${UTILS_TIMING_SRC_DIR}
    ${UTILS_TIMING_BIN_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/src/share/cxx
    )

  # Unfortunately, due to F90 dependencies, we end up including most of the f90
  # srcs, and one of them happens to include config.h explicitly
  # We already set all the config options by hand, so we can simply touch an empty config.h
  EXECUTE_PROCESS (COMMAND touch config.h WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  IF (USE_NUM_PROCS)
    SET (NUM_CPUS ${USE_NUM_PROCS})
  ELSE()
    SET (NUM_CPUS 1)
  ENDIF()
  cxx_unit_test (forcing_ut "${FORCING_UT_F90_SRCS}" "${FORCING_UT_CXX_SRCS}" "${FORCING_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
  TARGET_LINK_LIBRARIES(forcing_ut thetal_kokkos_ut_lib)

  # ### Caar functor unit test

  SET (CAAR_UT_CXX_SRCS
    ${THETA_UT_DIR}/caar_ut.cpp
    )

  SET (CAAR_UT_F90_SRCS
    ${THETA_UT_DIR}/caar_interface.F90
    ${THETA_UT_DIR}/thetal_test_interface.F90
    ${SHARE_UT_DIR}/geometry_interface.F90
    )

  SET (CAAR_UT_INCLUDE_DIRS
    ${SRC_THETA_DIR}/cxx
    ${SRC_SHARE_DIR}
    ${SRC_SHARE_DIR}/cxx
    ${THETA_UT_DIR}
    ${THETA_LIB_MODULE_DIR}
    ${UTILS_TIMING_SRC_DIR}
    ${UTILS_TIMING_BIN_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/src/share/cxx
    )

  IF (USE_NUM_PROCS)
    SET (NUM_CPUS ${USE_NUM_PROCS})
  ELSE()
    SET (NUM_CPUS 1)
  ENDIF()
  cxx_unit_test (caar_ut "${CAAR_UT_F90_SRCS}" "${CAAR_UT_CXX_SRCS}" "${CAAR_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
  TARGET_LINK_LIBRARIES(caar_ut thetal_kokkos_ut_lib)

  # ### Remap functor unit test

  SET (REMAP_THETA_UT_CXX_SRCS
    ${THETA_UT_DIR}/remap_theta_ut.cpp
    )

  SET (REMAP_THETA_UT_F90_SRCS
    ${THETA_UT_DIR}/remap_theta_interface.F90
    ${THETA_UT_DIR}/thetal_test_interface.F90
    ${SHARE_UT_DIR}/geometry_interface.F90
    )

  SET (REMAP_THETA_UT_INCLUDE_DIRS
    ${SRC_THETA_DIR}/cxx
    ${SRC_SHARE_DIR}
    ${SRC_SHARE_DIR}/cxx
    ${THETA_UT_DIR}
    ${THETA_LIB_MODULE_DIR}
    ${UTILS_TIMING_SRC_DIR}
    ${UTILS_TIMING_BIN_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/src/share/cxx
    )

  IF (USE_NUM_PROCS)
    SET (NUM_CPUS ${USE_NUM_PROCS})
  ELSE()
    SET (NUM_CPUS 1)
  ENDIF()
  cxx_unit_test (remap_theta_ut "${REMAP_THETA_UT_F90_SRCS}" "${REMAP_THETA_UT_CXX_SRCS}" "${REMAP_THETA_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
  TARGET_LINK_LIBRARIES(remap_theta_ut thetal_kokkos_ut_lib)
ENDIF ()

# ### DIRK functor unit test

SET (DIRK_UT_CXX_SRCS
  ${THETA_UT_DIR}/dirk_ut.cpp
)

SET (DIRK_UT_F90_SRCS
  ${THETA_UT_DIR}/dirk_interface.F90
  ${THETA_UT_DIR}/thetal_test_interface.F90
  ${SHARE_UT_DIR}/geometry_interface.F90
)

SET (DIRK_UT_INCLUDE_DIRS
  ${SRC_THETA_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${THETA_UT_DIR}
  ${THETA_LIB_MODULE_DIR}
  ${UTILS_TIMING_SRC_DIR}
  ${UTILS_TIMING_BIN_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (dirk_ut "${DIRK_UT_F90_SRCS}" "${DIRK_UT_CXX_SRCS}" "${DIRK_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
TARGET_LINK_LIBRARIES(dirk_ut thetal_kokkos_ut_lib)

# ### Compose semi-Lagrangian transport unit tests

SET (COMPOSE_UT_CXX_SRCS
  ${THETA_UT_DIR}/compose_ut.cpp
)

SET (COMPOSE_UT_F90_SRCS
  ${THETA_UT_DIR}/compose_interface.F90
  ${THETA_UT_DIR}/thetal_test_interface.F90
  ${SHARE_UT_DIR}/geometry_interface.F90
)

SET (COMPOSE_UT_INCLUDE_DIRS
  ${SRC_THETA_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/compose
  ${SRC_SHARE_DIR}/cxx
  ${THETA_UT_DIR}
  ${THETA_LIB_MODULE_DIR}
  ${UTILS_TIMING_SRC_DIR}
  ${UTILS_TIMING_BIN_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

SET (NUM_CPUS 1)
cxx_unit_test (compose_ut "${COMPOSE_UT_F90_SRCS}" "${COMPOSE_UT_CXX_SRCS}" "${COMPOSE_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
TARGET_LINK_LIBRARIES(compose_ut thetal_kokkos_ut_lib)

# ### GllFvRemap unit tests

SET (GLLFVREMAP_UT_CXX_SRCS
  ${THETA_UT_DIR}/gllfvremap_ut.cpp
)

SET (GLLFVREMAP_UT_F90_SRCS
  ${THETA_UT_DIR}/gllfvremap_interface.F90
  ${THETA_UT_DIR}/compose_interface.F90
  ${THETA_UT_DIR}/thetal_test_interface.F90
  ${SHARE_UT_DIR}/geometry_interface.F90
)

SET (GLLFVREMAP_UT_INCLUDE_DIRS
  ${SRC_THETA_DIR}/cxx
  ${SRC_SHARE_DIR}
  ${SRC_SHARE_DIR}/cxx
  ${THETA_UT_DIR}
  ${THETA_LIB_MODULE_DIR}
  ${UTILS_TIMING_SRC_DIR}
  ${UTILS_TIMING_BIN_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/src/share/cxx
)

IF (USE_NUM_PROCS)
  SET (NUM_CPUS ${USE_NUM_PROCS})
ELSE()
  SET (NUM_CPUS 1)
ENDIF()
cxx_unit_test (gllfvremap_ut "${GLLFVREMAP_UT_F90_SRCS}" "${GLLFVREMAP_UT_CXX_SRCS}" "${GLLFVREMAP_UT_INCLUDE_DIRS}" "${CONFIG_DEFINES}" ${NUM_CPUS})
TARGET_LINK_LIBRARIES(gllfvremap_ut thetal_kokkos_ut_lib)
cxx_unit_test_add_test(gllfvremap_planar_ut gllfvremap_ut ${NUM_CPUS} "hommexx -planar")
