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
  CreateUnitTest(spa_tests "${SPA_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP spa_tests_ut_np1_omp1)
endif()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/map_ne2np4_to_ne4np4_mono.nc
               ${CMAKE_CURRENT_BINARY_DIR}/map_ne2np4_to_ne4np4_mono.nc COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/spa_data_for_testing.nc
               ${CMAKE_CURRENT_BINARY_DIR}/spa_data_for_testing.nc COPYONLY)
