include(CMakeParseArguments) # Needed for backwards compatibility
include(EkatCreateUnitTest)
include(EkatUtils)

# Create a list containing a range of integers
function (CreateRange resultVar BEG END)
  set(options SKIP_FIRST SKIP_LAST)
  set(arg1v INC)
  set(argMv)
  cmake_parse_arguments(CR "${options}" "${arg1v}" "${argMv}" ${ARGN})

  # Compute beg/end/inc based on input args
  if (CR_SKIP_FIRST)
    math(EXPR BEG "${BEG}+1")
  endif()
  if (CR_SKIP_LAST)
    math(EXPR END "${END}-1")
  endif()
  if (NOT CR_INC)
    set (CR_INC 1)
  endif()

  # Sanity check
  if (NOT CR_INC GREATER 0)
    message (FATAL_ERROR "INC must be a positive integer")
  endif()
  if (BEG GREATER END)
    message (FATAL_ERROR "BEG is larger than END")
  endif()

  # Create range list
  set (res_list)
  set (N ${BEG})
  while (NOT N GREATER END)
    list (APPEND res_list ${N})
    math (EXPR N "${N}+${CR_INC}")
  endwhile()

  # Set in parent scope
  set (${resultVar} ${res_list} PARENT_SCOPE)
endfunction()

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

macro(SetVarDependingOnTestSize var_name val_short val_medium val_long)
  string (TOUPPER ${SCREAM_TEST_SIZE} size)
  if (size STREQUAL "SHORT")
    set(${var_name} ${val_short})
  elseif(size STREQUAL "MEDIUM")
    set(${var_name} ${val_medium})
  elseif(size STREQUAL "LONG")
    set(${var_name} ${val_long})
  else()
    message("Error! Unrecognized testing size '${size}'.")
    message("  Valid test size options: ${SCREAM_TEST_VALID_SIZES}")
    message(FATAL_ERROR "Aborting.")
  endif()
endmacro()

set(SCREAM_CUT_EXEC_OPTIONS ${CUT_EXEC_OPTIONS})
set(SCREAM_CUT_EXEC_1V_ARGS ${CUT_EXEC_1V_ARGS})
set(SCREAM_CUT_EXEC_MV_ARGS ${CUT_EXEC_MV_ARGS})

set(SCREAM_CUT_TEST_OPTIONS ${CUT_TEST_OPTIONS})
set(SCREAM_CUT_TEST_1V_ARGS ${CUT_TEST_1V_ARGS};MINIMUM_TEST_LEVEL)
set(SCREAM_CUT_TEST_MV_ARGS ${CUT_TEST_MV_ARGS})

#
# Adjust Ekat CUT options based on SCREAM specifics
#

# Scream always excludes the ekat test session since it has its own
list(REMOVE_ITEM SCREAM_CUT_EXEC_OPTIONS USER_DEFINED_TEST_SESSION)

###############################################################################
function(CreateUnitTestExec exec_name)
###############################################################################
  # Call Ekat function, with a couple of extra params
  EkatCreateUnitTestExec("${exec_name}" ${ARGN}
    USER_DEFINED_TEST_SESSION LIBS eamxx_test_support)
endfunction(CreateUnitTestExec)

###############################################################################
function(CreateUnitTestFromExec test_name test_exec)
###############################################################################
  cmake_parse_arguments(cutfe "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}" ${ARGN})
  CheckMacroArgs(CreateUnitTestFromExec cutfe "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}")

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

  if (cutfe_MINIMUM_TEST_LEVEL)
    if (cutfe_MINIMUM_TEST_LEVEL GREATER SCREAM_TEST_LEVEL)
      message("Test ${test_name} skipped due to insufficient test level")
      return()
    endif()
    unset(cutfe_MINIMUM_TEST_LEVEL) # Ekat does not take this option
  endif()

  separate_cut_arguments(cutfe "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}" options)

  EkatCreateUnitTestFromExec("${test_name}" "${test_exec}" ${options})

endfunction(CreateUnitTestFromExec)

###############################################################################
function(CreateUnitTest test_name)
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
  CreateUnitTestExec("${test_name}" ${options_ExecPhase})

  #------------------------------#
  #      Create Tests Phase      #
  #------------------------------#

  separate_cut_arguments(cut "${SCREAM_CUT_TEST_OPTIONS}" "${SCREAM_CUT_TEST_1V_ARGS}" "${SCREAM_CUT_TEST_MV_ARGS}" options_TestPhase)
  CreateUnitTestFromExec("${test_name}" "${test_name}" ${options_TestPhase})

endfunction(CreateUnitTest)

###############################################################################
function(GetInputFile src_path)
###############################################################################
  # Input files are NOT downloaded at configure time (a network hiccup would
  # leave us without a build tree). Instead, we register a custom target that
  # fetches the file at build time:
  #  - get_input_<file>: one target per file (created once, globally)
  #  - get_input_dir_<dir>: one target per calling directory, part of ALL, that
  #    depends on the files registered in that directory. This way, running
  #    'make' inside a test folder fetches the data that folder needs.
  #  - eamxx_input_data: umbrella target fetching every registered file.
  # We rely on CIME to download for us if this is a CIME case.
  if (SCREAM_CIME_BUILD)
    return()
  endif()
  if (NOT SCREAM_INPUT_ROOT)
    message(FATAL_ERROR "No SCREAM_INPUT_ROOT set, cannot try to locate input file ${src_path}")
  endif()

  string(MAKE_C_IDENTIFIER "${src_path}" file_id)
  set(file_tgt get_input_${file_id})
  if (NOT TARGET ${file_tgt})
    add_custom_target(${file_tgt}
      COMMAND ${CMAKE_COMMAND}
        -DINPUT_ROOT=${SCREAM_INPUT_ROOT}
        -DREL_PATH=${src_path}
        -DCHECK_INPUT_SCRIPT=${SCREAM_SRC_DIR}/../scripts/check-input
        -P ${SCREAM_SRC_DIR}/../cmake/DownloadInputFile.cmake
      COMMENT "Checking input file ${src_path}"
      VERBATIM)
    add_dependencies(eamxx_input_data ${file_tgt})
  endif()

  file(RELATIVE_PATH rel_dir ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  string(MAKE_C_IDENTIFIER "${rel_dir}" dir_id)
  set(dir_tgt get_input_dir_${dir_id})
  if (NOT TARGET ${dir_tgt})
    add_custom_target(${dir_tgt} ALL)
  endif()
  add_dependencies(${dir_tgt} ${file_tgt})
endfunction(GetInputFile)
