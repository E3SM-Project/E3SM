include(CMakeParseArguments) # Needed for backwards compatibility
include(EkatCreateUnitTest)

# This function takes the following arguments:
#    - target_name: the name of the executable
#    - target_srcs: a list of src files for the executable
#      Note: no need to include catch_main; this macro will add it
#    - scream_libs: a list of scream libraries needed by the executable
#    - config_defs (optional): a list of additional defines for the compiler, default is nothing
#    - mpi_ranks (optional): the number of mpi ranks for the test, if 2 values, it's a range, if 3, it's a range plus increment. default is np=1
#    - threads (optional): the number of threads for the test, if 2 values, it's a range, if 3, it's a range plus an increment. default is 1 thread
#      Note: One test will be created per combination of valid mpi-rank and thread value

function(CreateUnitTest target_name target_srcs scream_libs)
  set(options OPTIONAL EXCLUDE_MAIN_CPP)
  set(oneValueArgs EXE_ARGS DEP)
  set(multiValueArgs MPI_RANKS THREADS CONFIG_DEFS INCLUDE_DIRS COMPILER_FLAGS LABELS)
  cmake_parse_arguments(CreateUnitTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set (TEST_INCLUDE_DIRS
       ${SCREAM_INCLUDE_DIRS}
       ${CMAKE_CURRENT_SOURCE_DIR}
       ${CMAKE_CURRENT_BINARY_DIR}
       ${SCREAM_F90_MODULES}
       ${CreateUnitTest_INCLUDE_DIRS}
  )

  set (test_libs ${scream_libs})
  list(APPEND test_libs ${SCREAM_TPL_LIBRARIES})

  if (CreateUnitTest_EXCLUDE_MAIN_CPP)
    EkatCreateUnitTest(${target_name} "${target_srcs}"
      DEP "${CreateUnitTest_DEP}"
      MPI_EXEC_NAME ${SCREAM_MPIRUN_EXE}
      MPI_NP_FLAG -np
      MPI_RANKS "${CreateUnitTest_MPI_RANKS}"
      THREADS "${CreateUnitTest_THREADS}"
      EXE_ARGS "${CreateUnitTest_EXE_ARGS}"
      MPI_EXTRA_ARGS "${SCREAM_MPI_EXTRA_ARGS}"
      COMPILER_DEFS "${CreateUnitTest_CONFIG_DEFS}"
      INCLUDE_DIRS "${TEST_INCLUDE_DIRS}"
      COMPILER_FLAGS "${CreateUnitTest_COMPILER_FLAGS}"
      LIBS "${test_libs}"
      LIBS_DIRS "${SCREAM_TPL_LIBRARY_DIRS}"
      LINKER_FLAGS "${SCREAM_LINK_FLAGS}"
      LABELS "${CreateUnitTest_LABELS}"
      EXCLUDE_MAIN_CPP
    )
  else ()
    EkatCreateUnitTest(${target_name} "${target_srcs}"
      DEP "${CreateUnitTest_DEP}"
      MPI_EXEC_NAME ${SCREAM_MPIRUN_EXE}
      MPI_NP_FLAG -np
      MPI_RANKS "${CreateUnitTest_MPI_RANKS}"
      THREADS "${CreateUnitTest_THREADS}"
      EXE_ARGS "${CreateUnitTest_EXE_ARGS}"
      MPI_EXTRA_ARGS "${SCREAM_MPI_EXTRA_ARGS}"
      COMPILER_DEFS "${CreateUnitTest_CONFIG_DEFS}"
      INCLUDE_DIRS "${TEST_INCLUDE_DIRS}"
      COMPILER_FLAGS "${CreateUnitTest_COMPILER_FLAGS}"
      LIBS "${test_libs}"
      LIBS_DIRS "${SCREAM_TPL_LIBRARY_DIRS}"
      LINKER_FLAGS "${SCREAM_LINK_FLAGS}"
      LABELS "${CreateUnitTest_LABELS}"
    )
  endif ()

endfunction(CreateUnitTest)
