cmake_minimum_required(VERSION 3.9)

set(CTEST_BUILD_NAME "scream_unit_tests${BUILD_NAME_MOD}")

get_filename_component(working_dir ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)
set(CTEST_SOURCE_DIRECTORY "${working_dir}")
set(CTEST_BINARY_DIRECTORY "${working_dir}/ctest-build")

if(NOT DEFINED dashboard_model)
  set(dashboard_model Experimental)
endif()
If(NOT DEFINED dashboard_track)
  set(dashboard_track E3SM_SCREAM)
endif()

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(CTEST_CONFIGURE_COMMAND "${CMAKE_COMMAND} ${CTEST_SOURCE_DIRECTORY}")

ctest_start(${dashboard_model} TRACK ${dashboard_track})

ctest_configure()
ctest_build(FLAGS "-j8")

if (NOT BUILD_ONLY)
  ctest_test(RETURN_VALUE TEST_RESULTS PARALLEL_LEVEL ${PARALLEL_LEVEL})
  if (NOT TEST_RESULTS EQUAL 0)
    set(TEST_FAILS True)
  endif()
endif()

if (NOT NO_SUBMIT)
  ctest_submit(RETRY_COUNT 10 RETRY_DELAY 60)
endif()

if (TEST_FAILS)
  message(FATAL_ERROR "Test had fails")
endif()