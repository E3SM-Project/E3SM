SET(UNITTESTER_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

include (EkatBuildCatch2)
# Run a unit test executable with NUM_CPUS ranks. Extra args are passed to the executable.
# The MPI launcher is set up in homme_setup_testing (cmake/HommeTesting.cmake).
function(cxx_unit_test_add_test test_name target_name NUM_CPUS)
  SET(EXTRA_ARGS ${ARGN})
  ADD_TEST(NAME ${test_name}
           COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NUM_CPUS} ${HOMME_MPIEXEC_PREFLAGS_LIST} $<TARGET_FILE:${target_name}> ${HOMME_MPIEXEC_POSTFLAGS_LIST} ${EXTRA_ARGS})
  SET_TESTS_PROPERTIES(${test_name} PROPERTIES LABELS "unit" PROCESSORS ${NUM_CPUS})
endfunction()

macro(cxx_unit_test target_name target_f90_srcs target_cxx_srcs include_dirs config_defines NUM_CPUS)
  ADD_EXECUTABLE(${target_name} ${UNITTESTER_DIR}/tester.cpp ${target_f90_srcs} ${target_cxx_srcs})
  # Build this exec as part of the test-execs target
  ADD_DEPENDENCIES(test-execs ${target_name})

  cxx_unit_test_add_test(${target_name}_test ${target_name} ${NUM_CPUS})

  TARGET_LINK_LIBRARIES(${target_name} timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  TARGET_LINK_LIBRARIES(${target_name} Kokkos::kokkos)
  IF (HOMME_USE_MKL)
    IF (MKL_TYPE STREQUAL "oneMKL")
      TARGET_LINK_LIBRARIES(${target_name} -qmkl)
    ELSEIF (MKL_TYPE STREQUAL "Intel MKL")
      TARGET_LINK_LIBRARIES(${target_name} -mkl)
    ENDIF ()
  ENDIF()

  # Link csm_share lib
  TARGET_LINK_LIBRARIES(${target_name} csm_share)

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
  ELSE ()
    ADD_DEFINITIONS(-DHOMME_WITHOUT_PIOLIBRARY)
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

  TARGET_LINK_LIBRARIES(${target_name} Catch2::Catch2)
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${CATCH_INCLUDE}")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${include_dirs}")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${SRC_SHARE};${SRC_BASE};${UNITTESTER_DIR}")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${PIO_INCLUDE_DIRS}")
  TARGET_INCLUDE_DIRECTORIES(${target_name} PUBLIC "${CMAKE_BINARY_DIR}/src")

endmacro(cxx_unit_test)
