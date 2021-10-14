INCLUDE (ScreamUtils)

SET (NEED_LIBS spa scream_io physics_share scream_share)
set(SPA_TESTS_SRCS
    spa_unit_tests.cpp
#    spa_read_remap_from_file_test.cpp
#    spa_init_pio_for_test.cpp
#    spa_finish_pio_for_test.cpp
    spa_read_data_from_file_test.cpp
    ) # SPA_TESTS_SRCS

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(spa_tests "${SPA_TESTS_SRCS}" "${NEED_LIBS}" 
    LABELS "spa"
    MPI_RANKS 1 ${SCREAM_TEST_MAX_RANKS}
    PROPERTIES RESOURCE_LOCK "spa_unit_tests"
  )
endif()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_data_for_testing.nc
               ${CMAKE_CURRENT_BINARY_DIR}/spa_data_for_testing.nc COPYONLY)
