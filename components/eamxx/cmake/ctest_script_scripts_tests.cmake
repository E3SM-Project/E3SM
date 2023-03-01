cmake_minimum_required(VERSION 3.9)

set(CTEST_BUILD_NAME "scripts_tests")

set(CTEST_SOURCE_DIRECTORY "${SCREAM_ROOT}/scripts/tests")
set(CTEST_BINARY_DIRECTORY "${BUILD_WORK_DIR}")

if(NOT DEFINED dashboard_model)
  set(dashboard_model Experimental)
endif()
if(NOT DEFINED dashboard_track)
  set(dashboard_track E3SM_SCREAM)
endif()

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

ctest_start(${dashboard_model} TRACK ${dashboard_track})

ctest_configure(RETURN_VALUE CONFIG_ERROR_CODE)

if (CONFIG_ERROR_CODE)
  set (TEST_FAILS TRUE)
else ()

  ctest_test(RETURN_VALUE TEST_ERROR_CODE)

  if (TEST_ERROR_CODE)
    set(TEST_FAILS TRUE)
  endif()
endif ()

if (NOT NO_SUBMIT)
  ctest_submit(RETRY_COUNT 10 RETRY_DELAY 60)
endif()

if (TEST_FAILS)
  message(FATAL_ERROR "Test had fails")
endif()
