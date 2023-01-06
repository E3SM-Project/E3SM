# macro for building PREQx primitive equations executable

macro(preqx_setup)


INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
LINK_DIRECTORIES    ( ${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS} )

SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/preqx)
SET(SRC_DIR           ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR     ${HOMME_SOURCE_DIR}/src/share)
SET(TRILINOS_ZOLTAN_DIR ${HOMME_SOURCE_DIR}/src/zoltan)
SET(TEST_SRC_DIR      ${HOMME_SOURCE_DIR}/src/test_src)

SET(UTILS_TIMING_DIR  ${HOMME_BINARY_DIR}/utils/cime/CIME/non_py/src/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS} ${UTILS_TIMING_DIR})

# Find F90 files in target directory
FILE(GLOB TARGET_F90  ${TARGET_DIR}/*.F90 ${SRC_SHARE_DIR}/*.F90 ${SRC_DIR}/preqx/share/*.F90 ${TEST_SRC_DIR}/*.F90)


SET(PREQX_SRCS_F90
  ${TARGET_F90}
  ${SRC_DIR}/checksum_mod.F90 
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/common_movie_mod.F90
  ${SRC_DIR}/interpolate_driver_mod.F90
  ${SRC_DIR}/interp_movie_mod.F90
  ${SRC_DIR}/interpolate_driver_mod.F90
  ${SRC_DIR}/netcdf_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90
  ${SRC_DIR}/prim_main.F90 
  ${SRC_DIR}/prim_movie_mod.F90
  ${SRC_DIR}/prim_restart_mod.F90
  ${SRC_DIR}/repro_sum_mod.F90
  ${SRC_DIR}/restart_io_mod.F90
  ${SRC_DIR}/surfaces_mod.F90
  ${SRC_DIR}/test_mod.F90
)

IF (${HOMME_USE_TRILINOS})
SET (PREQX_SRCS_ZOLTAN   ${TRILINOS_ZOLTAN_DIR}/zoltan_interface.c  ${TRILINOS_ZOLTAN_DIR}/zoltan_cppinterface.cpp)
SET_SOURCE_FILES_PROPERTIES( ${TRILINOS_ZOLTAN_DIR}/zoltan_cppinterface.cpp PROPERTIES LANGUAGE CXX )
ENDIF()

IF (HOMMEXX_BFB_TESTING)
  # We need a couple more files/includes to correctly hook-up C routines
  SET(PREQX_SRCS_F90
      ${PREQX_SRCS_F90}
      ${SRC_SHARE_DIR}/cxx/utilities/bfb_mod.F90
  )

  SET(PREQX_SRCS_C
      ${PREQX_SRCS_C}
      ${SRC_SHARE_DIR}/cxx/utilities/BfbUtils.cpp
      ${SRC_SHARE_DIR}/cxx/ExecSpaceDefs.cpp
  )

  SET (EXEC_INCLUDE_DIRS
      ${EXEC_INCLUDE_DIRS}
      ${SRC_SHARE_DIR}/cxx
      ${HOMME_BINARY_DIR}/src/share/cxx/
  )
ENDIF()

# If the user specified a file for custom compiler options use those
IF (DEFINED PREQX_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(PREQX_CUSTOM_FLAGS_FILE PREQX_SRCS_F90)
ENDIF ()

# Make SRCS global so the tests can access it
SET(EXEC_SOURCES ${PREQX_SRCS_C} ${PREQX_SRCS_F90} ${PREQX_SRCS_ZOLTAN})

# Set up defaults
IF (NOT PREQX_NP)
  SET (PREQX_NP 4) 
ENDIF ()
IF (NOT PREQX_PLEV)
  SET (PREQX_PLEV 20) 
ENDIF ()
IF (NOT PREQX_USE_PIO)
  SET (PREQX_USE_PIO FALSE) 
ENDIF ()
IF (NOT PREQX_USE_ENERGY)
  SET (PREQX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PREQX_NC)
  SET (PREQX_NC 4)
ENDIF ()
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()
# for some of matt's old scripts which use preqx_qsize_d
IF (PREQX_QSIZE_D)      
  SET (QSIZE_D ${PREQX_QSIZE_D})
ENDIF ()

SET(USE_OPENACC FALSE)
IF (${ENABLE_HORIZ_OPENMP} AND ${ENABLE_COLUMN_OPENMP})
  SET(NESTED_OPENMP TRUE CACHE BOOL "Both horizontal and column threading")
ENDIF ()


endmacro(preqx_setup)

preqx_setup()
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(preqx preqx ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} ${PREQX_USE_PIO}  ${PREQX_USE_ENERGY} ${QSIZE_D})


