THETAL_KOKKOS_SETUP()

# This will be used to determine that we need to link to kokkos
SET(USE_KOKKOS_KERNELS ON)

# Set the variables for this test executable        NP  NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(theta-l-nlev72-kokkos theta-l_kokkos  4  4  72   FALSE   TRUE        10)

# Setting HOMME_TESTS_* variables, so the namelist.nl file in the exec 
# directory is usable. Since that namelist should be used for development
# and/or debugging purposes only, we make the test 'small' (ne=2, ndays=1),
# and pick qsize 4, rsplit 3 and moisture='notdry'.
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
SET (HOMME_TEST_NE 8)
SET (HOMME_TEST_TIME_STEP 450)
SET (HOMME_TEST_NDAYS 1)
SET (HOMME_TEST_QSIZE 10)
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
