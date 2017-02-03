# macro for building pese primitive equations executable

macro(pese_setup)

SET(TARGET_DIR        ${HOMME_SOURCE_DIR}/src/pese)
SET(UTILS_SHARE_DIR   ${HOMME_SOURCE_DIR}/utils/csm_share)
SET(SRC_DIR           ${HOMME_SOURCE_DIR}/src)
SET(SRC_SHARE_DIR     ${HOMME_SOURCE_DIR}/src/share)
SET(PHYSICS_DIR       ${HOMME_SOURCE_DIR}/physics)
SET(UTILS_TIMING_DIR  ${CMAKE_BINARY_DIR}/utils/cime/share/timing)

# Make INCLUDE_DIRS global so the tests can access it
SET (EXEC_INCLUDE_DIRS ${PIO_INCLUDE_DIRS} ${UTILS_TIMING_DIR} )

# Find F90 files in target directory
FILE(GLOB TARGET_F90  ${TARGET_DIR}/*.F90)
FILE(GLOB TEST_SRC_F90   ${SRC_SHARE_DIR}/*F90 ${SRC_SHARE_DIR}/vertical_se/*F90)

SET(PESE_SRCS_F90
  ${TARGET_F90}
  ${TEST_SRC_F90}
  ${SRC_DIR}/test_src/asp_tests.F90
  ${SRC_DIR}/test_src/baroclinic_inst_mod.F90
  ${SRC_DIR}/test_src/dcmip2012_test1_2_3.F90
  ${SRC_DIR}/test_src/dcmip2012_test4.F90
  ${SRC_DIR}/test_src/dcmip2012_test5.F90
  ${SRC_DIR}/test_src/held_suarez_mod.F90 
  ${SRC_DIR}/checksum_mod.F90
  ${SRC_DIR}/common_io_mod.F90
  ${SRC_DIR}/common_movie_mod.F90
  ${SRC_DIR}/interp_movie_mod.F90
  ${SRC_DIR}/interpolate_driver_mod.F90
  ${SRC_DIR}/netcdf_interp_mod.F90
  ${SRC_DIR}/netcdf_io_mod.F90
  ${SRC_DIR}/pio_io_mod.F90
  ${SRC_DIR}/prim_main.F90 
  ${SRC_DIR}/prim_movie_mod.F90
  ${SRC_DIR}/prim_restart_mod.F90
  ${SRC_DIR}/repro_sum_mod.F90
  ${SRC_DIR}/restart_io_mod.F90
  ${SRC_DIR}/surfaces_mod.F90
  ${SRC_DIR}/test_mod.F90
  ${UTILS_SHARE_DIR}/shr_const_mod.F90 
  ${UTILS_SHARE_DIR}/shr_file_mod.F90 
  ${UTILS_SHARE_DIR}/shr_kind_mod.F90
  ${UTILS_SHARE_DIR}/shr_mpi_mod.F90 
  ${UTILS_SHARE_DIR}/shr_spfn_mod.F90
  ${UTILS_SHARE_DIR}/shr_sys_mod.F90 
  ${UTILS_SHARE_DIR}/shr_vmath_mod.F90
)

SET(PESE_SRCS_C
  ${UTILS_SHARE_DIR}/shr_vmath_fwrap.c
)

# If the user specified a file for custom compiler options use those
IF (DEFINED PESE_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(PESE_CUSTOM_FLAGS_FILE PESE_SRCS_F90)
ENDIF ()

# Make SRCS global so the tests can access it
SET(EXEC_SOURCES ${PESE_SRCS} ${PESE_SRCS_C} ${PESE_SRCS_F90} )

# Set up defaults
IF (NOT PESE_NP)
  SET (PESE_NP 4) 
ENDIF ()
IF (NOT PESE_PLEV)
  SET (PESE_PLEV 20) 
ENDIF ()
IF (NOT PESE_USE_PIO)
  SET (PESE_USE_PIO FALSE) 
ENDIF ()
IF (NOT PESE_USE_ENERGY)
  SET (PESE_USE_ENERGY FALSE) 
ENDIF ()
IF (NOT PESE_NC)
  SET (PESE_NC 4)
ENDIF ()
IF (NOT QSIZE_D)
  SET (QSIZE_D 4)
ENDIF ()

SET(USE_OPENACC FALSE)
IF (${ENABLE_HORIZ_OPENMP} AND ${ENABLE_COLUMN_OPENMP})
  SET(NESTED_OPENMP TRUE CACHE BOOL "Both horizontal and column threading")
ENDIF ()

endmacro(pese_setup)

pese_setup()
############################################################################
# createTestExec(exec_name exec_type NP PLEV USE_PIO USE_ENERGY)
############################################################################
createTestExec(pese pese ${PESE_NP} ${PESE_NC} ${PESE_PLEV} ${PESE_USE_PIO}  ${PESE_USE_ENERGY} ${QSIZE_D})


