# Set the qsize_d in dimensions_mod.F90
SET(QSIZE_D 25)
SET(USE_OPENACC TRUE)
# Set the variables for this test executable
#                           NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(baroCam-acc preqx_acc 4   4   26   TRUE        TRUE)
