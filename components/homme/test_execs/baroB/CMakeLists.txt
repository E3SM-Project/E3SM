
# Set the variables for this test executable
#                          NP  NC PLEV USE_PIO WITH_ENERGY
createTestExec(baroB preqx  8   4   20    TRUE      FALSE)
IF (${HOMME_USE_TRILINOS})
INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
TARGET_LINK_LIBRARIES(baroB ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
ENDIF()

