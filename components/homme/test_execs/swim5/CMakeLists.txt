swim_setup()
#                          NP  NC PLEV USE_PIO  WITH_ENERGY  QSIZE_D
createTestExec(swim5 swim   4   4   1   FALSE        FALSE   4)

# we are doing somethign wrong here.  createTestExec should set these
# and they should not have to be repeated for every single test executable

TARGET_LINK_LIBRARIES(swim5 ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
SET_PROPERTY(TARGET swim5 PROPERTY LINKER_LANGUAGE Fortran)
