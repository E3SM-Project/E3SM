#_______________________________________________________________________
# create executable needed by DCMIP2016_test3 (interpolated output)

thetal_setup()
#              name                 target NP NC PLEV USE_PIO WITH_ENERGY QSIZE_D
createTestExec(theta-l-nlev40         theta-l  4  4  72   FALSE   TRUE       6)

# Setting HOMME_TESTS_* variables, so the namelist.nl file in the exec 
# directory is usable. Since that namelist should be used for development
# and/or debugging purposes only, we make the test 'small'
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
SET (HOMME_TEST_NE 8)
SET (HOMME_TEST_TIME_STEP 450)
SET (HOMME_TEST_NDAYS 1)

SET (HOMME_TEST_QSIZE 4)
SET (HOMME_TEST_RSPLIT 2)
SET (HOMME_TEST_NU 7e15 )
SET (HOMME_TEST_NUDIV 1e15)
SET (HOMME_TEST_NUTOP 2.5e5)
SET (HOMME_TEST_LIM 9)
SET (HOMME_TEST_HVSCALING 2)
SET (HOMME_TEST_HVS_TOM 0)
SET (HOMME_THETA_HY_MODE "true")
SET (HOMME_THETA_FORM 1)
SET (HOMME_TTYPE 5)

####note this needs renaming later
#SET (HOMME_TEST_HYDROSTATIC_MODE "true")

# Copy the needed input files to the binary dir
CONFIGURE_FILE (${CMAKE_SOURCE_DIR}/test/reg_test/namelists/theta.nl
                ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/movies)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/test/vcoord/${HOMME_TEST_VCOORD_INT_FILE}
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/test/vcoord/${HOMME_TEST_VCOORD_MID_FILE}
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)




