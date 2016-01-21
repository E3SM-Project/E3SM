
# Set the variables for this test executable
#                         NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(baroC preqx 4   4   26   FALSE        TRUE)
IF (${HOMME_USE_TRILINOS})
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES(baroC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
ENDIF()

