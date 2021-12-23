include(CMakeParseArguments) # Needed for backwards compatibility
include(EkatCreateUnitTest)
include(EkatUtils)

# This function takes the following arguments:
#    - test_name: the base name of the test. We create an executable with this name
#    - test_srcs: a list of src files for the executable.
#      Note: no need to include catch_main; this macro will add it
#    - scream_libs: a list of scream libraries needed by the executable
#    - compiler_defs (optional): a list of additional defines for the compiler, default is nothing
#    - compiler_flags (optional): a list of additional options for the compiler, default is nothing
#    - mpi_ranks (optional): the number of mpi ranks for the test, if 2 values, it's a range, if 3, it's a range plus increment. default is np=1
#    - threads (optional): the number of threads for the test, if 2 values, it's a range, if 3, it's a range plus an increment. default is 1 thread
# Notes:
#  - One test will be created per combination of valid mpi-rank and thread value
#  - compiler defs/flags can also be providedd on a per-language basis via COMPILER_[C|CXX|F]_[FLAGS|DEFS]

macro(SetVarDependingOnTestProfile var_name val_short val_medium val_long)
  string (TOUPPER ${SCREAM_TEST_PROFILE} profile)
  if (profile STREQUAL "SHORT")
    set(${var_name} ${val_short})
  elseif(profile STREQUAL "MEDIUM")
    set(${var_name} ${val_medium})
  elseif(profile STREQUAL "LONG")
    set(${var_name} ${val_long})
  else()
    message("Error! Unrecognized testing profile '${profile}'.")
    message("  Valid test profile options: ${SCREAM_TEST_VALID_PROFILES}")
    message(FATAL_ERROR "Aborting.")
  endif()
endmacro()

set(SCREAM_CUT_EXEC_OPTIONS ${CUT_EXEC_OPTIONS})
set(SCREAM_CUT_EXEC_1V_ARGS ${CUT_EXEC_1V_ARGS})
set(SCREAM_CUT_EXEC_MV_ARGS ${CUT_EXEC_MV_ARGS})

set(SCREAM_CUT_TEST_OPTIONS ${CUT_TEST_OPTIONS})
set(SCREAM_CUT_TEST_1V_ARGS ${CUT_TEST_1V_ARGS})
set(SCREAM_CUT_TEST_MV_ARGS ${CUT_TEST_MV_ARGS})

#
# Adjust Ekat CUT options based on SCREAM specifics
#

# Scream always excludes the ekat test session since it has its own
list(REMOVE_ITEM SCREAM_CUT_EXEC_OPTIONS EXCLUDE_TEST_SESSION)

# Libs are a position arg for SCREAM, not an optional arg like in EKAT
list(REMOVE_ITEM SCREAM_CUT_EXEC_MV_ARGS LIBS)

# MPI stuff is set in scream's cache config and is not configurable per-test
list(REMOVE_ITEM SCREAM_CUT_TEST_1V_ARGS MPI_EXEC_NAME MPI_NP_FLAG)
list(REMOVE_ITEM SCREAM_CUT_TEST_MV_ARGS MPI_EXTRA_ARGS)

###############################################################################
function(CreateUnitTestExec exec_name test_srcs scream_libs)
###############################################################################
  cmake_parse_arguments(cute "${SCREAM_CUT_EXEC_OPTIONS}" "${SCREAM_CUT_EXEC_1V_ARGS}" "${SCREAM_CUT_EXEC_MV_ARGS}" ${ARGN})
  CheckMacroArgs(CreateUnitTestExec cute "${SCREAM_CUT_EXEC_OPTIONS}" "${SCREAM_CUT_EXEC_1V_ARGS}" "${SCREAM_CUT_EXEC_MV_ARGS}")

  separate_cut_arguments(cute "${SCREAM_CUT_EXEC_OPTIONS}" "${SCREAM_CUT_EXEC_1V_ARGS}" "${SCREAM_CUT_EXEC_MV_ARGS}" options)

  set(TEST_INCLUDE_DIRS
    ${SCREAM_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${SCREAM_F90_MODULES}
  )

  set(test_libs "${scream_libs}")
  list(APPEND test_libs "${SCREAM_TPL_LIBRARIES}")

  if (SCREAM_Fortran_FLAGS)
    list(APPEND options COMPILER_F_FLAGS ${SCREAM_Fortran_FLAGS})
  endif ()

  EkatCreateUnitTestExec("${exec_name}" "${test_srcs}" ${options}
    EXCLUDE_TEST_SESSION LIBS ${test_libs} INCLUDE_DIRS ${TEST_INCLUDE_DIRS})

endfunction(CreateUnitTestExec)

###############################################################################
function(CreateUnitTestFromExec test_name test_exec)
###############################################################################
  cmake_parse_arguments(cutfe "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}" ${ARGN})
  CheckMacroArgs(CreateUnitTestExec cutfe "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}")

  separate_cut_arguments(cutfe "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}" options)

  if (SCREAM_MPI_EXTRA_ARGS)
    list(APPEND options MPI_EXTRA_ARGS ${SCREAM_MPI_EXTRA_ARGS})
  endif ()

  #
  # If asking for mpi/omp ranks/threads, verify we stay below the max number of threads
  #
  if (cutfe_MPI_RANKS OR cutfe_THREADS)
    list(LENGTH cutfe_MPI_RANKS NUM_MPI_RANK_ARGS)
    list(LENGTH cutfe_THREADS   NUM_THREAD_ARGS)
    if (NUM_MPI_RANK_ARGS EQUAL 0)
      set(MAX_RANKS 1)
    elseif (NUM_MPI_RANK_ARGS GREATER 1)
      list(GET cutfe_MPI_RANKS 1 MAX_RANKS)
    else()
      list(GET cutfe_MPI_RANKS 0 MAX_RANKS)
    endif()
    if (NUM_THREAD_ARGS EQUAL 0)
      set(MAX_THREADS 1)
    elseif (NUM_THREAD_ARGS GREATER 1)
      list(GET cutfe_THREADS   1 MAX_THREADS)
    else()
      list(GET cutfe_THREADS   0 MAX_THREADS)
    endif()
    math(EXPR NUM_THREADS_NEEDED ${MAX_RANKS}*${MAX_THREADS})
    if (${NUM_THREADS_NEEDED} GREATER ${SCREAM_TEST_MAX_TOTAL_THREADS})
      string (CONCAT msg
        "**************************************************************************\n"
        "Error! Invalid max threads/ranks combination. When invoking CreateUnitTest,\n"
        "you must ensure that the product of the max requested mpi ranks and\n"
        "the max requested omp threads is lower than SCREAM_TEST_MAX_TOTAL_THREADS.\n"
        "   - test exec name: ${test_name}\n"
        "   - requested max MPI ranks: ${MAX_RANKS}\n"
        "   - requested max OMP threads: ${MAX_THREADS}\n"
        "   - resulting max threads needed: ${NUM_THREADS_NEEDED}\n"
        "   - SCREAM_TEST_MAX_TOTAL_THREADS: ${SCREAM_TEST_MAX_TOTAL_THREADS}\n")
      message("${msg}")
      message(FATAL_ERROR "Aborting")
    endif()
  endif()

  EkatCreateUnitTestFromExec("${test_name}" "${test_exec}" ${options}
    MPI_EXEC_NAME ${SCREAM_MPIRUN_EXE} MPI_NP_FLAG ${SCREAM_MPI_NP_FLAG})

endfunction(CreateUnitTestFromExec)

###############################################################################
function(CreateUnitTest test_name test_srcs scream_libs)
###############################################################################
  set(options ${SCREAM_CUT_EXEC_OPTIONS} ${SCREAM_CUT_TEST_OPTIONS})
  set(oneValueArgs ${SCREAM_CUT_EXEC_1V_ARGS} ${SCREAM_CUT_TEST_1V_ARGS})
  set(multiValueArgs ${SCREAM_CUT_EXEC_MV_ARGS} ${SCREAM_CUT_TEST_MV_ARGS})

  # ecut = Ekat Create Unit Test
  cmake_parse_arguments(cut "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  CheckMacroArgs(CreateUnitTest cut "${options}" "${oneValueArgs}" "${multiValueArgs}")

  #------------------------------#
  #      Create Exec Phase       #
  #------------------------------#

  separate_cut_arguments(cut "${SCREAM_CUT_EXEC_OPTIONS}" "${SCREAM_CUT_EXEC_1V_ARGS}" "${SCREAM_CUT_EXEC_MV_ARGS}" options_ExecPhase)
  CreateUnitTestExec("${test_name}" "${test_srcs}" "${scream_libs}" ${options_ExecPhase})

  #------------------------------#
  #      Create Tests Phase      #
  #------------------------------#

  separate_cut_arguments(cut "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}" options_TestPhase)
  CreateUnitTestFromExec("${test_name}" "${test_name}" ${options_TestPhase})

endfunction(CreateUnitTest)
