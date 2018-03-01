# This function adds a pfunit test executable for BeTR.
function(add_betr_pfunit_test exe)

  # separate the input args into test sources and include dirs
  foreach (arg ${ARGN})
    if (arg MATCHES ".pfunit")
      list(APPEND test_sources ${arg})
    else()
      list(APPEND include_dirs ${arg})
    endif()
  endforeach()

  # call pfunit preprocessor to generate F90 sources
  include(PreprocessMacro)
  add_pfunit_sources(sources ${test_sources})

  include_directories(${CMAKE_BINARY_DIR}/mod)
  # add user specified include dirs necessary to find mod files for
  # compiling tests.
  foreach (dir ${include_dirs})
    include_directories(${CMAKE_BINARY_DIR}/${dir})
  endforeach()

  # need a copy of the pfunit driver
  configure_file(${CMAKE_BINARY_DIR}/include/driver.F90 driver.F90 COPYONLY)

  # create the executable
  add_executable(${exe} ${sources} driver.F90)
  target_link_libraries(${exe} pfunit;${BETR_LIBRARIES})
  set_target_properties(${exe} PROPERTIES LINKER_LANGUAGE Fortran)
  set_target_properties(${exe} PROPERTIES COMPILE_FLAGS "-DCMAKE_CURRENT_SOURCE_DIR=\\\"${CMAKE_CURRENT_SOURCE_DIR}\\\"")

  # define the test
  add_test(${exe} ${exe})
  set_tests_properties(${exe} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endfunction()
