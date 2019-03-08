set(SCREAM_TEST_SRCS
  ${SCREAM_SRC_DIR}/share/util/scream_catch_main.cpp
  scream_tests.cpp
)

# Create the dynamics lib
SET (PREQX_NP 4)
SET (PREQX_PLEV 72)
SET (PREQX_QSIZE_D 40)
SET (PREQX_USE_PIO FALSE)
SET (PREQX_USE_ENERGY FALSE)
SET (HOMME_TARGET "preqx")
CreateDynamicsLib()

link_directories(${SCREAM_LIBRARY_DIRS} ${SCREAM_TPL_LIBRARY_DIRS})

add_executable(scream_tests ${SCREAM_TEST_SRCS})
target_include_directories(scream_tests PUBLIC ${SCREAM_INCLUDE_DIRS} ${SCREAM_TPL_INCLUDE_DIRS} ${CATCH_INCLUDE_DIR})
target_link_libraries(scream_tests scream_control scream_dynamics p3 rrtmgp shoc scream_share ${SCREAM_TPL_LIBRARIES})
set_target_properties(scream_tests PROPERTIES LINK_FLAGS "${SCREAM_LINK_FLAGS}")

# Setting HOMME_TESTS_* variables, so the namelist.nl file in the exec 
# directory is usable. Since that namelist should be used for development
# and/or debugging purposes only, we make the test 'small' (ne=2, ndays=1),
# and pick qsize 4, rsplit 3 and moisture='notdry'.
SET (HOMME_TEST_VCOORD_INT_FILE acme-72i.ascii)
SET (HOMME_TEST_VCOORD_MID_FILE acme-72m.ascii)
SET (HOMME_TEST_NE 2)
SET (HOMME_TEST_NDAYS 1)
SET (HOMME_TEST_QSIZE 4)
SET (HOMME_TEST_RSPLIT 3)
SET (HOMME_TEST_MOISTURE notdry)

# Copy the needed input files to the binary dir
CONFIGURE_FILE (${HOMME_SOURCE_DIR}/test/reg_test/namelists/preqx.nl
                ${CMAKE_CURRENT_BINARY_DIR}/namelist.nl)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/movies)

FILE (MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/vcoord)

CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/vcoord/${HOMME_TEST_VCOORD_INT_FILE}
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)
CONFIGURE_FILE(${HOMME_SOURCE_DIR}/test/vcoord/${HOMME_TEST_VCOORD_MID_FILE}
               ${CMAKE_CURRENT_BINARY_DIR}/vcoord COPYONLY)

add_test(run_scream_tests scream_tests)
