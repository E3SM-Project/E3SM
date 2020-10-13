INCLUDE (ScreamUtils)

# Get or create the dynamics lib
#                 HOMME_TARGET   NP PLEV QSIZE_D USE_PIO USE_ENERGY 
CreateDynamicsLib("preqx_kokkos"  4   72   35    FALSE   FALSE)

SET (EXTRA_CONFIG_DEFS
     HAVE_CONFIG_H
     HOMMEXX_CONFIG_IS_CMAKE
)
SET (p3Lib p3)

SET (NEED_LIBS ${p3Lib} ${dynLibName} scream_control scream_share pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
SET (UNIT_TEST_EXTRA_LIBRARY_DIRS ${dynLibDir})

# Test atmosphere processes
CreateUnitTest(homme_p3 "homme_p3.cpp" "${NEED_LIBS}"
    COMPILER_DEFS "${EXTRA_CONFIG_DEFS}"
    EXE_ARGS "< namelist.nl"
    LABELS "dynamics;p3;driver")

# Set homme's test options, so that we can configure the namelist correctly
SET (HOMME_TEST_NE 2)
SET (HOMME_TEST_NDAYS 1)
SET (HOMME_TEST_QSIZE 4)
SET (HOMME_TEST_RSPLIT 3)
SET (HOMME_TEST_TIME_STEP 300)
SET (HOMME_TEST_MOISTURE notdry)
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)


# Configure the namelist into the test directory
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/reg_test/namelists/preqx.nl
               ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

# Create movies and vcoord directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/movies)
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

# Copy vcoord files to build dir
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/vcoord/acme-72i.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/vcoord/acme-72m.ascii
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)

# Copy p3 lookup tables to local data directory
FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_lookup_table_2.dat-v4
               ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE(${SCREAM_DATA_DIR}/p3_universal_constants.inp
               ${CMAKE_CURRENT_BINARY_DIR}/data COPYONLY)
