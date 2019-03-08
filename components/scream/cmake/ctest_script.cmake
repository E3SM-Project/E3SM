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

ctest_configure(CAPTURE_CMAKE_ERROR configure_ctest_error RETURN_VALUE configure_error)
message("configure error       = ${configure_error}")
message("configure ctest error = ${configure_ctest_error}")

ctest_build(CAPTURE_CMAKE_ERROR build_ctest_error RETURN_VALUE build_error FLAGS "-j8")
message("build error           = ${build_error}")
message("build ctest error     = ${build_ctest_error}")

if (NOT BUILD_ONLY)
  ctest_test(PARALLEL_LEVEL ${PARALLEL_LEVEL} CAPTURE_CMAKE_ERROR test_ctest_error RETURN_VALUE test_error)
  message("test error            = ${test_error}")
  message("test ctest error      = ${test_ctest_error}")
endif()

if (NOT NO_SUBMIT)
  ctest_submit(RETRY_COUNT 10 RETRY_DELAY 60 CAPTURE_CMAKE_ERROR submit_ctest_error RETURN_VALUE submit_error)
  message("submit error            = ${submit_error}")
  message("submit ctest error      = ${submit_ctest_error}")
endif()
