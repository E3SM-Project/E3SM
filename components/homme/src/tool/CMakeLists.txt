# macro for building the tools executable

macro(tool_setup)

SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/tool)
SET(UTILS_SHARE_DIR   ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR           ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR     ${HOMME_SOURCE_DIR}/src/share)
SET(TEST_SRC_DIR      ${HOMME_SOURCE_DIR}/src/test_src)
SET(UTILS_TIMING_DIR  ${HOMME_BINARY_DIR}/utils/cime/CIME/non_py/src/timing)
# Use theta-l to define the dycore-specific data structures.
SET(DYCORE_DIR        ${HOMME_SOURCE_DIR}/src/theta-l)
SET(DYCORE_SHARE_DIR  ${HOMME_SOURCE_DIR}/src/theta-l/share)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS} ${UTILS_TIMING_DIR} )

# Find F90 files in share and test directories
FILE(GLOB TARGET_F90 ${TARGET_DIR}/*.F90 ${DYCORE_SHARE_DIR}/*.F90 ${DYCORE_DIR}/*.F90 ${SRC_SHARE_DIR}/*.F90 ${TEST_SRC_DIR}/*.F90)

SET(TOOL_SRCS_F90
  ${TARGET_F90}
  ${SRC_DIR}/checksum_mod.F90 
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/common_movie_mod.F90
  ${SRC_DIR}/interpolate_driver_mod.F90
  ${SRC_DIR}/interp_movie_mod.F90
  ${SRC_DIR}/netcdf_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90
  ${SRC_DIR}/prim_movie_mod.F90
  ${SRC_DIR}/theta_restart_mod.F90
  ${SRC_DIR}/repro_sum_mod.F90
  ${SRC_DIR}/restart_io_mod.F90
  ${SRC_DIR}/surfaces_mod.F90
  ${SRC_DIR}/test_mod.F90
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90 
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90 
  ${UTILS_SHARE_DIR}/shr_file_mod.F90 
  ${UTILS_SHARE_DIR}/shr_vmath_mod.F90 
  ${UTILS_SHARE_DIR}/shr_const_mod.F90 
  ${UTILS_SHARE_DIR}/shr_spfn_mod.F90
)

# Make SRCS global so the tests can access it
SET(EXEC_SOURCES ${TOOL_SRCS_F90})

# Set up defaults
IF (NOT PREQX_NP)
  SET (PREQX_NP 4) 
ENDIF ()
IF (NOT PREQX_PLEV)
  SET (PREQX_PLEV 20) 
ENDIF ()
#IF (NOT PREQX_USE_PIO)
  SET (PREQX_USE_PIO TRUE)
#ENDIF ()
IF (NOT PREQX_USE_ENERGY)
  SET (PREQX_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PREQX_NC)
  SET (PREQX_NC 0)
ENDIF ()
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()
endmacro(tool_setup)

tool_setup()
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(homme_tool tool ${PREQX_NP} ${PREQX_NC} ${PREQX_PLEV} ${PREQX_USE_PIO} ${PREQX_USE_ENERGY} ${QSIZE_D})

