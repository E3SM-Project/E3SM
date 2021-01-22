SET(UNITTESTER_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

macro(cxx_unit_test target_name target_f90_srcs target_cxx_srcs include_dirs config_defines NUM_CPUS)
  ADD_EXECUTABLE(${target_name} ${UNITTESTER_DIR}/tester.cpp ${target_f90_srcs} ${target_cxx_srcs})
  #add exec to test_execs ,baseline, and check targets in makefile
  ADD_DEPENDENCIES(test-execs ${target_name})
  ADD_DEPENDENCIES(baseline ${target_name})
  ADD_DEPENDENCIES(check ${target_name})
  #IF (${NUM_CPUS} EQUAL 1)
  #  ADD_TEST(${target_name}_test ${target_name})
  #ELSE()
    ADD_TEST(${target_name}_test ${USE_MPIEXEC} -n ${NUM_CPUS} ${MPI_OPTIONS} "./${target_name}")
  #ENDIF()
  SET_TESTS_PROPERTIES(${target_name}_test PROPERTIES LABELS "unit")
  IF(BUILD_HOMME_WITHOUT_PIOLIBRARY)
    TARGET_COMPILE_DEFINITIONS(${target_name} PUBLIC HOMME_WITHOUT_PIOLIBRARY)
  ENDIF()

  TARGET_LINK_LIBRARIES(${target_name} timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} kokkos)

  STRING(TOUPPER "${PERFORMANCE_PROFILE}" PERF_PROF_UPPER)
  IF ("${PERF_PROF_UPPER}" STREQUAL "VTUNE")
    TARGET_LINK_LIBRARIES(${target_name} ittnotify)
  ENDIF()

  IF(UNIX AND NOT APPLE)
    TARGET_LINK_LIBRARIES(${target_name} rt)
  ENDIF()
  IF(NOT BUILD_HOMME_WITHOUT_PIOLIBRARY)
    IF(HOMME_USE_SCORPIO)
      TARGET_LINK_LIBRARIES(${target_name} piof pioc)
    ELSE ()
      TARGET_LINK_LIBRARIES(${target_name} pio)
    ENDIF ()
  ENDIF ()

  IF (NOT "${target_f90_srcs}" EQUAL "")
    # Fortran modules
    SET(MODULE_DIR ${CMAKE_CURRENT_BINARY_DIR}/${target_name}_modules)
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR})
  ENDIF()

  # Compile definitions and flags
  SET_TARGET_PROPERTIES(${target_name} PROPERTIES COMPILE_DEFINITIONS "${config_defines}")
  IF(BUILD_HOMME_WITHOUT_PIOLIBRARY)
    TARGET_COMPILE_DEFINITIONS(${target_name} PUBLIC HOMME_WITHOUT_PIOLIBRARY)
  ENDIF()

  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${HOMME_SOURCE_DIR}/test/unit_tests/catch2/include")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${include_dirs}")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${SRC_SHARE};${SRC_BASE};${UNITTESTER_DIR}")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${PIO_INCLUDE_DIRS};${UTILS_TIMING_DIR}")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${CMAKE_BINARY_DIR}/src")

endmacro(cxx_unit_test)
