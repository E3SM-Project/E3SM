PREQX_KOKKOS_SETUP()

# This will be used to determine that we need to link to kokkos
SET(USE_KOKKOS_KERNELS ON)

# Set the variables for this test executable
#                          NP  NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(preqx-nlev26-kokkos preqx_kokkos 4 4 26 FALSE FALSE 4)

# Setting HOMME_TESTS_* variables, so the namelist.nl file in the exec 
# directory is usable. Since that namelist should be used for development
# and/or debugging purposes only, we make the test 'small' (ne=2, ndays=1),
# and pick qsize 4, rsplit 3 and moisture='notdry'.
SET (HOMME_TEST_VCOORD_INT_FILE cami-26.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE camm-26.ascii)
SET (HOMME_TEST_NE 2)
SET (HOMME_TEST_NDAYS 1)
SET (HOMME_TEST_QSIZE 4)
SET (HOMME_TEST_RSPLIT 3)
SET (HOMME_TEST_MOISTURE notdry)
# Copy the needed input files to the binary dir
CONFIGURE_FILE (${CMAKE_SOURCE_DIR}/test/reg_test/namelists/preqx.nl
                ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/movies)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/test/vcoord/camm-26.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/test/vcoord/cami-26.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
