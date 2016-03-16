# Set the qsize_d in dimensions_mod.F90
SET(QSIZE_D 25)
SET(USE_OPENACC TRUE)
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OPENACC_Fortran_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENACC_Linker_FLAGS}")
# If the user specified a file for custom compiler options use those
IF (DEFINED PREQX_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(PREQX_CUSTOM_FLAGS_FILE PREQX_SRCS_F90)
ENDIF ()
# Set the variables for this test executable
#                           NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(baroCam-acc preqx_acc 4   4   26   TRUE        TRUE)
