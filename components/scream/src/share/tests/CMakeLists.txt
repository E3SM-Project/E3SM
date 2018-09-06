link_directories(${SCREAM_LIBRARY_DIRS})
add_executable(util_unittest_f90 unit_test.f90)
add_executable(util_tests
  util_tests.cpp
  ../util/scream_catch_main.cpp)
foreach (exe util_tests util_unittest_f90)
  target_link_libraries(${exe} scream_share ${SCREAM_LIBRARIES})
  target_include_directories(${exe} PUBLIC ${SCREAM_INCLUDE_DIRS} ${CATCH_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
  set_target_properties(${exe} PROPERTIES
    LINK_FLAGS "${SCREAM_LINK_FLAGS}"
    Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
  add_test(run_${exe} ${exe})
endforeach ()


