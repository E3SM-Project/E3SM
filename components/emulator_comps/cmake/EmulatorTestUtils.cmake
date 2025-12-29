#==============================================================================
# Emulator Test Utilities
#
# CMake macros for creating and running tests.
#==============================================================================

include(CMakeParseArguments)

#------------------------------------------------------------------------------
# Test Levels
#------------------------------------------------------------------------------
set(EMULATOR_TEST_LEVEL_UNIT 1)
set(EMULATOR_TEST_LEVEL_COMPONENT 2)
set(EMULATOR_TEST_LEVEL_SYSTEM 3)

if(NOT DEFINED EMULATOR_TEST_LEVEL)
  set(EMULATOR_TEST_LEVEL 3 CACHE STRING "Test level (1=unit, 2=component, 3=system)")
endif()

message(STATUS "Test level: ${EMULATOR_TEST_LEVEL} (1=unit, 2=component, 3=system)")

#------------------------------------------------------------------------------
# EmulatorCreateTest
#
# Creates a test executable and registers it with CTest.
#
# Usage:
#   EmulatorCreateTest(test_name
#     SOURCES source1.cpp source2.cpp
#     LIBS lib1 lib2
#     LABELS label1 label2
#     MPI_RANKS 1 2 4
#     LEVEL 1|2|3
#   )
#------------------------------------------------------------------------------
function(EmulatorCreateTest test_name)
  cmake_parse_arguments(ECT "" "LEVEL" "SOURCES;LIBS;LABELS;MPI_RANKS" ${ARGN})
  
  # Check test level - skip if requested level is higher than configured
  if(ECT_LEVEL AND ECT_LEVEL GREATER EMULATOR_TEST_LEVEL)
    message(STATUS "  Skipping ${test_name} (level ${ECT_LEVEL} > ${EMULATOR_TEST_LEVEL})")
    return()
  endif()
  
  # Validate sources
  if(NOT ECT_SOURCES)
    message(FATAL_ERROR "EmulatorCreateTest: SOURCES is required for ${test_name}")
  endif()
  
  # Create executable
  add_executable(${test_name} ${ECT_SOURCES})
  
  # Include Catch2
  if(DEFINED EMULATOR_CATCH2_INCLUDE_DIR)
    target_include_directories(${test_name} PRIVATE ${EMULATOR_CATCH2_INCLUDE_DIR})
  endif()
  
  # Link test support library and any additional libs
  target_link_libraries(${test_name} PRIVATE
    emulator_test_support
    ${ECT_LIBS}
  )
  
  target_compile_features(${test_name} PRIVATE cxx_std_17)
  
  # Build labels list
  set(test_labels ${ECT_LABELS})
  if(ECT_LEVEL EQUAL 1)
    list(APPEND test_labels "unit")
  elseif(ECT_LEVEL EQUAL 2)
    list(APPEND test_labels "component")
  elseif(ECT_LEVEL EQUAL 3)
    list(APPEND test_labels "system")
  endif()
  
  # Register tests with CTest
  if(ECT_MPI_RANKS)
    # Option to skip MPI tests (useful for login nodes where srun isn't available)
    if(EMULATOR_SKIP_MPI_TESTS)
      message(STATUS "  Skipping MPI test: ${test_name} (EMULATOR_SKIP_MPI_TESTS=ON)")
      # Still add a serial version
      add_test(
        NAME ${test_name}
        COMMAND ${test_name}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )
      set_tests_properties(${test_name} PROPERTIES 
        LABELS "${test_labels}"
        ENVIRONMENT "OMP_NUM_THREADS=1"
      )
    else()
      # Create test for each MPI rank count
      foreach(nranks ${ECT_MPI_RANKS})
        set(full_test_name "${test_name}_np${nranks}")
        add_test(
          NAME ${full_test_name}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nranks} 
                  $<TARGET_FILE:${test_name}>
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        set_tests_properties(${full_test_name} PROPERTIES
          LABELS "${test_labels};mpi"
          PROCESSORS ${nranks}
          ENVIRONMENT "OMP_NUM_THREADS=1"
        )
      endforeach()
    endif()
  else()
    # Serial test
    add_test(
      NAME ${test_name} 
      COMMAND ${test_name}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(${test_name} PROPERTIES 
      LABELS "${test_labels}"
      ENVIRONMENT "OMP_NUM_THREADS=1"
    )
  endif()
  
  message(STATUS "  Added test: ${test_name} (level ${ECT_LEVEL})")
  
endfunction()

#------------------------------------------------------------------------------
# Convenience macros for each test level
#------------------------------------------------------------------------------

# Unit test (level 1) - Fast, isolated function/class tests
macro(EmulatorUnitTest test_name)
  EmulatorCreateTest(${test_name} LEVEL 1 ${ARGN})
endmacro()

# Component test (level 2) - Single component with minimal dependencies
macro(EmulatorComponentTest test_name)
  EmulatorCreateTest(${test_name} LEVEL 2 ${ARGN})
endmacro()

# System test (level 3) - Full system integration tests
macro(EmulatorSystemTest test_name)
  EmulatorCreateTest(${test_name} LEVEL 3 ${ARGN})
endmacro()
