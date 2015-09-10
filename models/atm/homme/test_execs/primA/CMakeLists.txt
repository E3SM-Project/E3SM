INCLUDE_DIRECTORIES ( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS} )
# Set the variables for this test executable
#                         NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(primA prim  4   4   26   FALSE        TRUE)
TARGET_LINK_LIBRARIES(primA ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
SET_PROPERTY(TARGET primA PROPERTY LINKER_LANGUAGE Fortran)
