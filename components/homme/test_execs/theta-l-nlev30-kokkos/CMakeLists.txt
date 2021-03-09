THETAL_KOKKOS_SETUP()

# Set the variables for this test executable        NP  NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(theta-l-nlev30-kokkos theta-l_kokkos  4  4  30   FALSE   TRUE        6)

# Setting HOMME_TESTS_* variables, so the namelist.nl file in the exec 
# directory is usable. Since that namelist should be used for development
# and/or debugging purposes only, we make the test 'small' (ne=2, ndays=1),
# and pick qsize 4, rsplit 3 and moisture='notdry'.
SET (HOMME_TEST_VCOORD_INT_FILE cami-30.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE camm-30.ascii)
SET (HOMME_TEST_NE 8)
SET (HOMME_TEST_TIME_STEP 450)
SET (HOMME_TEST_NDAYS 1)
SET (HOMME_TEST_QSIZE 4)
SET (HOMME_TEST_RSPLIT 2)
SET (HOMME_TEST_QSPLIT 2)
SET (HOMME_TEST_MOISTURE dry)
SET (HOMME_TEST_HYDROSTATIC_MODE "true")

# Copy the needed input files to the binary dir
CONFIGURE_FILE (${CMAKE_SOURCE_DIR}/test/reg_test/namelists/thetanh.nl
                ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/movies)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/test/vcoord/${HOMME_TEST_VCOORD_INT_FILE}
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/test/vcoord/${HOMME_TEST_VCOORD_MID_FILE}
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
