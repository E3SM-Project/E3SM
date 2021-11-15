INCLUDE (ScreamUtils)

SET (NEED_LIBS spa scream_io physics_share scream_share)

CreateUnitTest(spa_read_data_test "spa_read_data_from_file_test.cpp" "${NEED_LIBS}"
  LABELS "spa"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)
CreateUnitTest(spa_read_remap_test "spa_read_remap_from_file_test.cpp" "${NEED_LIBS}"
  LABELS "spa"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)
CreateUnitTest(spa_one_to_one_remap_test "spa_one_to_one_remap_test.cpp" "${NEED_LIBS}"
  LABELS "spa"
  MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_data_for_testing.nc
               ${CMAKE_CURRENT_BINARY_DIR}/spa_data_for_testing.nc COPYONLY)
