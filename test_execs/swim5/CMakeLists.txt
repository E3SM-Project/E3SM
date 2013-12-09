INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
# Set the variables for this test executable
#                          NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(swim5 swim   4   4   1   FALSE        FALSE)
TARGET_LINK_LIBRARIES(swim5 ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})

