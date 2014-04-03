# Set the qsize_d in dimensions_mod.F90
SET(QSIZE_D 25)
# Set the variables for this test executable
#                            NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(perfTest preqx 4   4   26   FALSE        TRUE)
