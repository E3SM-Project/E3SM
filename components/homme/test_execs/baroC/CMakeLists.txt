
# If the user specified a file for custom compiler options use those
IF (DEFINED PREQX_CUSTOM_FLAGS_FILE)
  setCustomCompilerFlags(PREQX_CUSTOM_FLAGS_FILE PREQX_SRCS_F90)
ENDIF ()
# Set the variables for this test executable
#                         NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(baroC preqx 4   4   26   FALSE        TRUE)
