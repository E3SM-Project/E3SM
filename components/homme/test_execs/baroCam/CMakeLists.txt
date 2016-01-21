
# Set the qsize_d in dimensions_mod.F90
SET(QSIZE_D 25)
# Set the variables for this test executable
#                           NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(baroCam preqx 4   4   26   TRUE        TRUE)

IF (${HOMME_USE_TRILINOS})
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES(baroCam ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
ENDIF()

