# Utilities to create Homme's ctest tests.
#
# Everything here is plain CMake: tests are registered with add_test. A system test 'foo' runs
#   <mpi-launcher> <exe> foo.nl foo.log
# in its folder ${HOMME_TESTS_DIR}/foo, which contains the input files and the outputs (foo.log,
# movies/*.nc). The checks on the outputs are separate tests, which run CMake scripts (in cmake/testing)
# via 'cmake -P'. Use 'ctest -R foo -V' to see the commands of foo and its checks.
#
# The functions below require homme_setup_testing() to have been called.

set (HOMME_TESTING_SCRIPTS_DIR ${CMAKE_CURRENT_LIST_DIR}/testing)

# Set up the variables used by all the tests (MPI launcher, baseline dirs, ...)
macro (homme_setup_testing)
  ##############################################
  #               MPI launcher                 #
  ##############################################
  # Tests are launched with
  #   ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} <nranks> [${HOMME_MPIEXEC_THREADS_FLAG} <nthreads>] ${MPIEXEC_PREFLAGS} <exe> ${MPIEXEC_POSTFLAGS}
  # The MPIEXEC_XYZ names are the standard ones used by CMake's FindMPI module, which also
  # tries to set some defaults for them. The threads flag is only used for runs where the test
  # specifies a number of threads (e.g., '-c' for srun, which is the default if the launcher is srun).
  # NOTE: the tests must run inside an allocation. There is no support for submitting jobs to a queue.
  if (USE_QUEUING)
    string (CONCAT msg
      "Error! Submitting Homme tests to a queue is no longer supported (USE_QUEUING=${USE_QUEUING}).\n"
      "  Run the tests from inside an allocation instead, and set the MPI launcher\n"
      "  via MPIEXEC_EXECUTABLE (e.g., srun) and MPIEXEC_PREFLAGS (e.g., '-K --cpu-bind=cores').")
    message (FATAL_ERROR "${msg}")
  endif()

  # Deprecated names, still accepted for user-provided machine files
  if (DEFINED USE_MPIEXEC AND NOT "${USE_MPIEXEC}" STREQUAL "")
    message (DEPRECATION "USE_MPIEXEC is deprecated, use MPIEXEC_EXECUTABLE instead.")
    set (MPIEXEC_EXECUTABLE "${USE_MPIEXEC}" CACHE FILEPATH "Executable for running MPI programs" FORCE)
  endif()
  if (DEFINED USE_MPI_OPTIONS AND NOT "${USE_MPI_OPTIONS}" STREQUAL "")
    message (DEPRECATION "USE_MPI_OPTIONS is deprecated, use MPIEXEC_PREFLAGS instead.")
    set (MPIEXEC_PREFLAGS "${USE_MPI_OPTIONS}" CACHE STRING "Flags for MPIEXEC_EXECUTABLE, before the executable" FORCE)
  endif()
  if (DEFINED USE_MPI_RUN_SCRIPT AND NOT "${USE_MPI_RUN_SCRIPT}" STREQUAL "")
    message (FATAL_ERROR "USE_MPI_RUN_SCRIPT is no longer supported. Use MPIEXEC_EXECUTABLE and MPIEXEC_PREFLAGS instead.")
  endif()

  if (NOT MPIEXEC_EXECUTABLE)
    # Not set by the user, and FindMPI (if called) could not find it either
    set (MPIEXEC_EXECUTABLE "mpiexec" CACHE FILEPATH "Executable for running MPI programs" FORCE)
  endif()
  set (MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "Flag used by MPIEXEC_EXECUTABLE to specify the number of ranks")
  set (MPIEXEC_PREFLAGS     "" CACHE STRING "Flags for MPIEXEC_EXECUTABLE, before the executable")
  set (MPIEXEC_POSTFLAGS    "" CACHE STRING "Flags for MPIEXEC_EXECUTABLE, after the executable")
  get_filename_component (_homme_mpiexec_name "${MPIEXEC_EXECUTABLE}" NAME)
  if (_homme_mpiexec_name STREQUAL "srun")
    set (_homme_threads_flag_default "-c")
  else()
    set (_homme_threads_flag_default "")
  endif()
  set (HOMME_MPIEXEC_THREADS_FLAG "${_homme_threads_flag_default}" CACHE STRING
       "Flag used by MPIEXEC_EXECUTABLE to specify the threads per rank, for tests that set the number of threads")
  # MPIEXEC_PREFLAGS/POSTFLAGS can be CMake lists (as in FindMPI, e.g. "--bind-to;core") or
  # space-separated strings (as in our machine files, e.g. "-K --cpu-bind=cores"), or a mix of the two
  homme_split_launcher_flags (HOMME_MPIEXEC_PREFLAGS_LIST "${MPIEXEC_PREFLAGS}")
  homme_split_launcher_flags (HOMME_MPIEXEC_POSTFLAGS_LIST "${MPIEXEC_POSTFLAGS}")
  message (STATUS "Homme tests MPI launcher: ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} <N> ${MPIEXEC_PREFLAGS} <exe> ${MPIEXEC_POSTFLAGS}")

  ##############################################
  #             Max ranks per test             #
  ##############################################
  if (USE_NUM_PROCS)
    # USE_NUM_PROCS is the number of ranks CIME has allocated (e.g. from a
    # _P96 test option). Treat it as the cap: each test uses as many ranks
    # as it needs, but no more than USE_NUM_PROCS.
    set (HOMME_TESTS_MAX_RANKS ${USE_NUM_PROCS})
    message (STATUS "Homme tests max ranks (from USE_NUM_PROCS): ${USE_NUM_PROCS}")
  else ()
    # The default of 24 matches a typical CIME _P24 test layout
    set (HOMME_TESTING_MAX_PROCS 24 CACHE STRING "Maximum MPI ranks per test when USE_NUM_PROCS is not set")
    set (HOMME_TESTS_MAX_RANKS ${HOMME_TESTING_MAX_PROCS})
    message (STATUS "Homme tests max ranks (HOMME_TESTING_MAX_PROCS): ${HOMME_TESTING_MAX_PROCS}")
  endif ()

  # Per-test time limit (seconds). 0 means no limit.
  set (HOMME_TESTING_TIMELIMIT 1800 CACHE STRING "Time limit (in seconds) of each Homme system test (0 = no limit)")

  ##############################################
  #                 Baselines                  #
  ##############################################
  # Tests run in ${HOMME_TESTS_DIR}/<test-name>
  set (HOMME_TESTS_DIR ${CMAKE_BINARY_DIR}/tests)

  # The baselines for test X are in ${HOMME_BASELINES_DIR}/X/movies.
  # NOTE: CIME's HOMME test relies on the '/tests/baseline' suffix.
  if (DEFINED HOMME_BASELINE_DIR)
    set (HOMME_BASELINES_DIR ${HOMME_BASELINE_DIR}/tests/baseline)
  else ()
    set (HOMME_BASELINES_DIR ${HOMME_TESTS_DIR}/baseline)
  endif ()
  message (STATUS "Homme baselines dir: ${HOMME_BASELINES_DIR}")

  # Allow to skip the comparison of the tests' output against baselines (e.g., if we only care about
  # comparisons within the same build, like f90 vs cxx)
  option (HOMME_SKIP_BASELINE_COMPARE "Skip the comparison against baseline results in the tests" OFF)

  # If ON, instead of tests comparing against baselines, create tests that copy the output files
  # in the baselines dir. These tests have the 'baseline_gen' label, so baselines can be generated
  # by running 'ctest -L baseline_gen'.
  option (HOMME_GENERATE_BASELINES "Create tests that copy the tests outputs to the baselines dir" OFF)
  if (HOMME_GENERATE_BASELINES AND HOMME_SKIP_BASELINE_COMPARE)
    message (FATAL_ERROR "HOMME_GENERATE_BASELINES and HOMME_SKIP_BASELINE_COMPARE cannot be both ON")
  endif()

  # The 'baseline' target runs the tests that generate baselines, and copies their outputs in
  # ${HOMME_TESTS_DIR}/baseline. This script, created at the end of the configuration, does the copy.
  set (HOMME_COPY_BASELINES_SCRIPT ${HOMME_TESTS_DIR}/copy-baselines.cmake)
  set_property (GLOBAL PROPERTY HOMME_COPY_BASELINES_CMDS "")
endmacro()

# Split launcher flags into a list of arguments: first split the CMake list, then split each
# entry as a unix command line (so that quotes can be used for arguments containing spaces).
function (homme_split_launcher_flags out_var flags)
  set (args)
  foreach (entry IN LISTS flags)
    separate_arguments (entry_args UNIX_COMMAND "${entry}")
    list (APPEND args ${entry_args})
  endforeach()
  set (${out_var} "${args}" PARENT_SCOPE)
endfunction()

# Create a system test, that is, a test that runs a Homme executable with one or more namelists
#
# Usage:
#   homme_add_system_test (<name>
#     EXEC <exe-target>          # The Homme executable to run (must be a target)
#     RANKS <n>                  # Number of MPI ranks (adjusted to the available ranks, see homme_resolve_test_ranks)
#     NAMELISTS <files...>       # Namelists, run in the given order. They are configured (via configure_file) in the test dir.
#     [MAX_RANKS <n>]            # Max number of ranks that the test supports
#     [VCOORD_FILES <files...>]  # Copied in <test-dir>/vcoord
#     [MESH_FILES <files...>]    # Copied in <test-dir>
#     [NCL_FILES <files...>]     # Copied in <test-dir>
#     [OUTPUTS <files...>]       # Files in <test-dir>/movies to compare against the baselines
#     [OMP_NAMELISTS <files...>] # If ENABLE_HORIZ_OPENMP=ON, run these after NAMELISTS, with threads (test <name>_omp)
#     [OMP_THREADS <n>]          # Number of threads per rank for OMP_NAMELISTS (the ranks are RANKS/OMP_THREADS)
#     [REF_OUTPUTS <files...>]   # Compare these files in <test-dir>/movies against the CHECK_OUTPUTS (test <name>_ref_cmp)
#     [CHECK_OUTPUTS <files...>]
#     [REF_TOL <tol>]            # Tolerance on the normalized RMS difference for the REF/CHECK comparison (default: 0)
#     [LABELS <labels...>]
#     [TIMEOUT <seconds>])       # Default: HOMME_TESTING_TIMELIMIT
#
# The test is not created if EXEC is not a target. The following tests are created:
#   <name>              : runs EXEC with the namelist (see homme_add_run_tests). With N>1 namelists, there are N tests
#                         <name>_run1, ..., <name>_runN, run in this order.
#                         Labels: system, the model of EXEC, and baseline_gen (if OUTPUTS are given).
#   <name>_baseline_cmp : compares the OUTPUTS against the baselines (if OUTPUTS are given). Label: baseline_cmp.
#                         Not created if HOMME_SKIP_BASELINE_COMPARE or HOMME_GENERATE_BASELINES is ON.
#   <name>_baseline_gen : copies the OUTPUTS in the baselines dir (if HOMME_GENERATE_BASELINES is ON). Label: baseline_gen
#   <name>_omp          : runs the OMP_NAMELISTS after the NAMELISTS (if given, and ENABLE_HORIZ_OPENMP is ON)
#   <name>_ref_cmp      : compares CHECK_OUTPUTS against REF_OUTPUTS (if given)
# The tests that check the outputs of <name> require the fixture <name> (set up by all the runs of
# NAMELISTS), so 'ctest -R <name>_baseline_cmp' also runs <name> (or <name>_run1, ..., <name>_runN) first.
function (homme_add_system_test name)
  set (options)
  set (oneValueArgs EXEC RANKS MAX_RANKS OMP_THREADS REF_TOL TIMEOUT)
  set (multiValueArgs NAMELISTS VCOORD_FILES MESH_FILES NCL_FILES OUTPUTS OMP_NAMELISTS REF_OUTPUTS CHECK_OUTPUTS LABELS)
  cmake_parse_arguments (PARSE_ARGV 1 hast "${options}" "${oneValueArgs}" "${multiValueArgs}")
  if (hast_UNPARSED_ARGUMENTS)
    message (FATAL_ERROR "homme_add_system_test(${name}): unrecognized arguments: ${hast_UNPARSED_ARGUMENTS}")
  endif()
  foreach (arg IN ITEMS EXEC RANKS NAMELISTS)
    if (NOT hast_${arg})
      message (FATAL_ERROR "homme_add_system_test(${name}): missing required argument ${arg}")
    endif()
  endforeach()

  if (NOT TARGET ${hast_EXEC})
    message (STATUS "Not configuring test ${name} since it depends upon the executable ${hast_EXEC}, which isn't built with this configuration")
    return()
  endif()
  message (STATUS "Adding test: ${name}, using exec ${hast_EXEC}")

  # Allow wildcards in file names. The namelists are run in the given order, so only expand
  # the entries that contain a wildcard (file(GLOB) sorts the files).
  foreach (arg IN ITEMS VCOORD_FILES MESH_FILES NCL_FILES)
    file (GLOB hast_${arg} ${hast_${arg}})
  endforeach()
  foreach (arg IN ITEMS NAMELISTS OMP_NAMELISTS)
    set (files)
    foreach (f IN LISTS hast_${arg})
      if (f MATCHES "[*?]")
        file (GLOB expanded ${f})
        list (APPEND files ${expanded})
      else()
        list (APPEND files ${f})
      endif()
    endforeach()
    set (hast_${arg} ${files})
  endforeach()

  #################################
  #     Set up the test folder    #
  #################################
  set (test_dir ${HOMME_TESTS_DIR}/${name})
  foreach (nl IN LISTS hast_NAMELISTS hast_OMP_NAMELISTS)
    # Some namelists contain cmake variables (e.g., to generate multiple testcases
    # with different ne or ndays), so use configure_file to replace them with their values
    get_filename_component (nl_name ${nl} NAME)
    configure_file (${nl} ${test_dir}/${nl_name})
  endforeach()
  if (hast_VCOORD_FILES)
    file (COPY ${hast_VCOORD_FILES} DESTINATION ${test_dir}/vcoord)
  endif()
  if (hast_MESH_FILES OR hast_NCL_FILES)
    file (COPY ${hast_MESH_FILES} ${hast_NCL_FILES} DESTINATION ${test_dir})
  endif()
  file (MAKE_DIRECTORY ${test_dir}/movies ${test_dir}/restart)

  #################################
  #        The run test(s)        #
  #################################
  set (TEST_NAME ${name}) # Used in homme_resolve_test_ranks messages
  homme_resolve_test_ranks (nranks ${hast_RANKS} "${hast_MAX_RANKS}" "${hast_NAMELISTS}")

  get_target_property (model ${hast_EXEC} HOMME_MODEL)
  set (labels system ${hast_LABELS})
  if (model)
    list (APPEND labels ${model})
  endif()
  if (hast_OUTPUTS)
    list (APPEND labels baseline_gen)
  endif()
  if (DEFINED hast_TIMEOUT)
    set (timeout ${hast_TIMEOUT})
  else()
    set (timeout ${HOMME_TESTING_TIMELIMIT})
  endif()

  homme_add_run_tests (${name} ${hast_EXEC} ${test_dir} ${nranks} 0 "${hast_NAMELISTS}" "" "${labels}" ${timeout})

  set (all_fixtures ${name})
  if (hast_OMP_NAMELISTS AND ENABLE_HORIZ_OPENMP)
    if (NOT hast_OMP_THREADS)
      message (FATAL_ERROR "homme_add_system_test(${name}): OMP_NAMELISTS requires OMP_THREADS")
    endif()
    math (EXPR omp_mod "${nranks} % ${hast_OMP_THREADS}")
    if (NOT omp_mod EQUAL 0)
      message (STATUS "In test ${name} NUM_CPUS not divisible by OMP_NUM_THREADS.")
    endif()
    math (EXPR omp_nranks "${nranks} / ${hast_OMP_THREADS}")
    # The OMP runs may need files from the first runs (e.g., restart files)
    homme_add_run_tests (${name}_omp ${hast_EXEC} ${test_dir} ${omp_nranks} ${hast_OMP_THREADS}
                         "${hast_OMP_NAMELISTS}" ${name} "system;${model};${hast_LABELS}" ${timeout})
    list (APPEND all_fixtures ${name}_omp)
  endif()

  #################################
  #   Checks on the test outputs  #
  #################################
  # NOTE: lists passed to add_test would be expanded into separate arguments, so we join them with commas
  string (REPLACE ";" "," outputs "${hast_OUTPUTS}")
  string (REPLACE ";" "," ref_outputs "${hast_REF_OUTPUTS}")
  string (REPLACE ";" "," check_outputs "${hast_CHECK_OUTPUTS}")
  if (hast_OUTPUTS)
    if (HOMME_GENERATE_BASELINES)
      add_test (NAME ${name}_baseline_gen
                COMMAND ${CMAKE_COMMAND} -DSRC_DIR=${test_dir}/movies -DDST_DIR=${HOMME_BASELINES_DIR}/${name}/movies
                        -DFILES=${outputs}
                        -P ${HOMME_TESTING_SCRIPTS_DIR}/HommeCopyFiles.cmake)
      set_tests_properties (${name}_baseline_gen PROPERTIES FIXTURES_REQUIRED ${name} LABELS "baseline_gen;${hast_LABELS}")
    elseif (NOT HOMME_SKIP_BASELINE_COMPARE)
      add_test (NAME ${name}_baseline_cmp
                COMMAND ${CMAKE_COMMAND} -DCPRNC=${CPRNC_BINARY} -DSRC_DIR=${test_dir}/movies
                        -DTGT_DIR=${HOMME_BASELINES_DIR}/${name}/movies -DFILES=${outputs}
                        -DLOG_PREFIX=${name}
                        -P ${HOMME_TESTING_SCRIPTS_DIR}/HommeCprncCompare.cmake
                WORKING_DIRECTORY ${test_dir})
      set_tests_properties (${name}_baseline_cmp PROPERTIES FIXTURES_REQUIRED ${name} LABELS "baseline_cmp;${hast_LABELS}")
    endif()

    # Used by the 'baseline' target
    set_property (GLOBAL APPEND PROPERTY HOMME_COPY_BASELINES_CMDS
      "homme_copy_files(\"${test_dir}/movies\" \"${HOMME_TESTS_DIR}/baseline/${name}/movies\" \"${outputs}\")")
  endif()

  if (hast_CHECK_OUTPUTS)
    list (LENGTH hast_REF_OUTPUTS nref)
    list (LENGTH hast_CHECK_OUTPUTS nchk)
    if (nchk GREATER nref)
      message (FATAL_ERROR "homme_add_system_test(${name}): CHECK_OUTPUTS has more entries than REF_OUTPUTS")
    endif()
    if (NOT DEFINED hast_REF_TOL)
      set (hast_REF_TOL 0.0)
    endif()
    add_test (NAME ${name}_ref_cmp
              COMMAND ${CMAKE_COMMAND} -DCPRNC=${CPRNC_BINARY} -DTGT_DIR=${test_dir}/movies -DSRC_DIR=${test_dir}/movies
                      -DTGT_FILES=${ref_outputs} -DFILES=${check_outputs} -DTOL=${hast_REF_TOL}
                      -DLOG_PREFIX=${name}.ref
                      -P ${HOMME_TESTING_SCRIPTS_DIR}/HommeCprncCompare.cmake
              WORKING_DIRECTORY ${test_dir})
    set_tests_properties (${name}_ref_cmp PROPERTIES FIXTURES_REQUIRED "${all_fixtures}" LABELS "system;${hast_LABELS}")
  endif()

  add_dependencies (test-execs ${hast_EXEC})
endfunction()

# Internal: create the tests that run an executable, one test per namelist. With one namelist,
# the test is <name>, otherwise the tests are <name>_run1, <name>_run2, ... Each test runs (in test_dir)
#   ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} <nranks> [${HOMME_MPIEXEC_THREADS_FLAG} <nthreads>]
#       ${MPIEXEC_PREFLAGS} <exe> ${MPIEXEC_POSTFLAGS} <namelist> <test>.log
# so Homme reads the namelist file, and writes its log in <test_dir>/<test>.log.
# The runs happen in order (each run requires the previous one, and the first one requires the
# fixture <requires>, if not empty), and they are all setup tests of the fixture <name>, so the
# tests requiring <name> run after all of them. If nthreads>0, OMP_NUM_THREADS is set for the runs.
function (homme_add_run_tests name exe test_dir nranks nthreads namelists requires labels timeout)
  set (cmd ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nranks})
  if (nthreads GREATER 0 AND HOMME_MPIEXEC_THREADS_FLAG)
    list (APPEND cmd ${HOMME_MPIEXEC_THREADS_FLAG} ${nthreads})
  endif()
  list (APPEND cmd ${HOMME_MPIEXEC_PREFLAGS_LIST} $<TARGET_FILE:${exe}> ${HOMME_MPIEXEC_POSTFLAGS_LIST})

  set (nprocs ${nranks})
  if (nthreads GREATER 0)
    math (EXPR nprocs "${nranks} * ${nthreads}")
  endif()

  list (LENGTH namelists nruns)
  set (irun 0)
  foreach (nl IN LISTS namelists)
    math (EXPR irun "${irun} + 1")
    if (nruns EQUAL 1)
      set (test_name ${name})
    else()
      set (test_name ${name}_run${irun})
    endif()
    get_filename_component (nl_name ${nl} NAME)
    add_test (NAME ${test_name} COMMAND ${cmd} ${nl_name} ${test_name}.log WORKING_DIRECTORY ${test_dir})
    set (fixtures ${name} ${test_name})
    list (REMOVE_DUPLICATES fixtures)
    set_tests_properties (${test_name} PROPERTIES FIXTURES_SETUP "${fixtures}"
                          LABELS "${labels}" PROCESSORS ${nprocs})
    if (requires)
      set_tests_properties (${test_name} PROPERTIES FIXTURES_REQUIRED ${requires})
    endif()
    if (timeout GREATER 0)
      set_tests_properties (${test_name} PROPERTIES TIMEOUT ${timeout})
    endif()
    if (nthreads GREATER 0)
      set_tests_properties (${test_name} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=${nthreads};OMP_STACKSIZE=128M")
    endif()
    set (requires ${test_name})
  endforeach()
endfunction()

# Compute the number of ranks for a test, given its requested number of ranks (nranks_in)
# and the max number of ranks it supports (max_ranks, can be empty):
#  - lower it to min(HOMME_TESTS_MAX_RANKS, max_ranks)
#  - raise it toward that cap if there is headroom: use the largest divisor of the number
#    of elements that is <= the cap (so that no rank gets 0 elements). The number of
#    elements is deduced from the namelists (ne, or ne_x/ne_y for planar grids).
function (homme_resolve_test_ranks out_var nranks_in max_ranks namelists)
  set (NUM_CPUS ${nranks_in})
  set (_HOMME_MAX_PROCS ${HOMME_TESTS_MAX_RANKS})
  if (max_ranks AND max_ranks LESS _HOMME_MAX_PROCS)
    set (_HOMME_MAX_PROCS ${max_ranks})
  endif ()
  if (NUM_CPUS GREATER _HOMME_MAX_PROCS)
    message (STATUS "Test ${TEST_NAME}: lowering NUM_CPUS from ${NUM_CPUS} to ${_HOMME_MAX_PROCS} (MAX_NUM_PROCS=${HOMME_TESTS_MAX_RANKS}, MAX_NUM_CPUS=${max_ranks})")
    set (NUM_CPUS ${_HOMME_MAX_PROCS})
  endif ()
  # Scan all namelists for ne/ne_x/ne_y, using the first file that contains a non-zero value
  # (a mesh-file namelist with "ne = 0" may sort ahead of the one with a non-zero ne).
  # Both hard-coded values (e.g. ne = 10) and ne = ${HOMME_TEST_NE} are handled.
  if (NUM_CPUS LESS _HOMME_MAX_PROCS AND namelists)
    set (_HOMME_TOTAL_ELEM 0)
    foreach (_HOMME_NL_FILE IN LISTS namelists)
      if (EXISTS "${_HOMME_NL_FILE}" AND _HOMME_TOTAL_ELEM EQUAL 0)
        file (READ "${_HOMME_NL_FILE}" _HOMME_NL_CONTENTS)
        if ("${_HOMME_NL_CONTENTS}" MATCHES "ne_x[ \t]*=[ \t]*([0-9]+)")
          # Planar
          set (_HOMME_NE_X "${CMAKE_MATCH_1}")
          if ("${_HOMME_NL_CONTENTS}" MATCHES "ne_y[ \t]*=[ \t]*([0-9]+)")
            set (_HOMME_NE_Y "${CMAKE_MATCH_1}")
          else ()
            # ne_y not found; assume square domain
            set (_HOMME_NE_Y "${_HOMME_NE_X}")
            message (STATUS "Test ${TEST_NAME}: ne_y not found in namelist; assuming ne_y = ne_x = ${_HOMME_NE_X} for MPI-count computation")
          endif ()
          math (EXPR _HOMME_TOTAL_ELEM "${_HOMME_NE_X} * ${_HOMME_NE_Y}")
        else ()
          # Cube-sphere. Skip ne = 0, which means the namelist uses a mesh file
          if ("${_HOMME_NL_CONTENTS}" MATCHES "[\n\t ]ne[ \t]*=[ \t]*([0-9]+)")
            set (_HOMME_NE_VAL "${CMAKE_MATCH_1}")
          elseif ("${_HOMME_NL_CONTENTS}" MATCHES "[\n\t ]ne[ \t]*=[ \t]*\\$\\{HOMME_TEST_NE\\}")
            set (_HOMME_NE_VAL "${HOMME_TEST_NE}")
          else ()
            set (_HOMME_NE_VAL "")
          endif ()
          if (_HOMME_NE_VAL)
            math (EXPR _HOMME_TOTAL_ELEM "6 * ${_HOMME_NE_VAL} * ${_HOMME_NE_VAL}")
          endif ()
        endif ()
      endif ()
    endforeach ()
    if (_HOMME_TOTAL_ELEM EQUAL 0)
      # Unknown grid size (e.g., RRM mesh file): keep the requested NUM_CPUS.
      message (STATUS "Test ${TEST_NAME}: keeping NUM_CPUS at ${NUM_CPUS} (total_elem unknown/mesh-file, max_procs=${_HOMME_MAX_PROCS})")
    elseif (_HOMME_TOTAL_ELEM GREATER_EQUAL _HOMME_MAX_PROCS)
      # Large grid: every rank gets >= 1 element
      message (STATUS "Test ${TEST_NAME}: raising NUM_CPUS from ${NUM_CPUS} to ${_HOMME_MAX_PROCS} (total_elem=${_HOMME_TOTAL_ELEM}, max_procs=${_HOMME_MAX_PROCS})")
      set (NUM_CPUS ${_HOMME_MAX_PROCS})
    else ()
      # Small grid: find the largest divisor of total_elem that is <= _HOMME_MAX_PROCS
      set (_HOMME_BEST_CPUS ${NUM_CPUS})
      set (_HOMME_CANDIDATE ${_HOMME_TOTAL_ELEM})
      while (_HOMME_CANDIDATE GREATER NUM_CPUS)
        math (EXPR _HOMME_REM "${_HOMME_TOTAL_ELEM} % ${_HOMME_CANDIDATE}")
        if (_HOMME_REM EQUAL 0)
          set (_HOMME_BEST_CPUS ${_HOMME_CANDIDATE})
          break ()
        endif ()
        math (EXPR _HOMME_CANDIDATE "${_HOMME_CANDIDATE} - 1")
      endwhile ()
      if (_HOMME_BEST_CPUS GREATER NUM_CPUS)
        message (STATUS "Test ${TEST_NAME}: raising NUM_CPUS from ${NUM_CPUS} to ${_HOMME_BEST_CPUS} (total_elem=${_HOMME_TOTAL_ELEM}, max_procs=${_HOMME_MAX_PROCS})")
        set (NUM_CPUS ${_HOMME_BEST_CPUS})
      endif ()
    endif ()
  endif ()
  set (${out_var} ${NUM_CPUS} PARENT_SCOPE)
endfunction()

# Create a test that compares the output of a F90 test and a C++ test, both bit-for-bit
# (netcdf files) and via the diagnostics written in the logs (only for theta tests).
#   homme_add_cxx_vs_f90_test (<name> F90_TEST <f90-test> CXX_TEST <cxx-test> OUTPUTS <files...> [LABELS ...])
function (homme_add_cxx_vs_f90_test name)
  cmake_parse_arguments (PARSE_ARGV 1 hacvf "" "F90_TEST;CXX_TEST" "OUTPUTS;LABELS")
  string (REPLACE ";" "," outputs "${hacvf_OUTPUTS}")
  add_test (NAME ${name}
            COMMAND ${CMAKE_COMMAND} -DCPRNC=${CPRNC_BINARY}
                    -DF90_DIR=${HOMME_TESTS_DIR}/${hacvf_F90_TEST} -DF90_TEST_NAME=${hacvf_F90_TEST}
                    -DCXX_DIR=${HOMME_TESTS_DIR}/${hacvf_CXX_TEST} -DCXX_TEST_NAME=${hacvf_CXX_TEST}
                    -DFILES=${outputs}
                    -P ${HOMME_TESTING_SCRIPTS_DIR}/HommeCxxVsF90.cmake
            WORKING_DIRECTORY ${HOMME_TESTS_DIR}/${hacvf_CXX_TEST})
  set_tests_properties (${name} PROPERTIES FIXTURES_REQUIRED "${hacvf_F90_TEST};${hacvf_CXX_TEST}"
                        LABELS "cxx_vs_f90;${hacvf_LABELS}")
endfunction()

# Check prescribed-wind tracer transport test error-norm output, which is available with
# some tests. See cmake/testing/HommeCheckErrorNorms.cmake for details.
#   homme_add_error_norms_test (<test-name> <anchor> <log-file> <upper-bounds>)
function (homme_add_error_norms_test test_name anchor filename upper_bounds)
  add_test (NAME "${test_name}_l2err"
            COMMAND ${CMAKE_COMMAND} -DANCHOR=${anchor} -DFILENAME=${filename} -DUPPER_BOUNDS=${upper_bounds}
                    -P ${HOMME_TESTING_SCRIPTS_DIR}/HommeCheckErrorNorms.cmake
            WORKING_DIRECTORY ${HOMME_TESTS_DIR}/${test_name})
  set_tests_properties ("${test_name}_l2err" PROPERTIES FIXTURES_REQUIRED "${test_name}" LABELS "system")
endfunction()

# Create the script used by the 'baseline' target to copy the tests outputs in ${HOMME_TESTS_DIR}/baseline
function (homme_write_copy_baselines_script)
  get_property (cmds GLOBAL PROPERTY HOMME_COPY_BASELINES_CMDS)
  set (content "# Generated by CMake: copies the outputs of the tests in ${HOMME_TESTS_DIR}/baseline\n")
  string (APPEND content "include (\"${HOMME_TESTING_SCRIPTS_DIR}/HommeCopyFiles.cmake\")\n")
  foreach (cmd IN LISTS cmds)
    string (APPEND content "${cmd}\n")
  endforeach()
  file (WRITE ${HOMME_COPY_BASELINES_SCRIPT} "${content}")
endfunction()

###############################################################################
#  Support for the test description files in test/reg_test/run_tests/*.cmake  #
###############################################################################

# Read a test description file (which sets TEST_NAME, EXEC_NAME, NUM_CPUS, NAMELIST_FILES, ...),
# and create the corresponding test(s). If PROFILE is set, the test is created for that
# testing profile (see set_homme_tests_parameters), and its name is changed accordingly.
macro (createTest testFile)
  set (THIS_TEST_INPUT ${HOMME_SOURCE_DIR}/test/reg_test/run_tests/${testFile})

  resetTestVariables()
  set (HOMME_ROOT ${HOMME_SOURCE_DIR})

  include (${THIS_TEST_INPUT})

  set (THIS_TEST_LABELS)
  if (DEFINED PROFILE)
    set_homme_tests_parameters(${testFile} ${PROFILE})
    if ("${TEST_NAME}" MATCHES "theta-f")
      set (TEST_NAME "${TEST_NAME}-ne${HOMME_TEST_NE}-nu${HOMME_TEST_NU}-ndays${HOMME_TEST_NDAYS}")
    elseif ("${TEST_NAME}" MATCHES "preqx-nlev")
      set (TEST_NAME "${TEST_NAME}-ne${HOMME_TEST_NE}-ndays${HOMME_TEST_NDAYS}")
    endif()
    set (THIS_TEST_LABELS ${PROFILE})
  endif ()

  if ("${NUM_CPUS}" STREQUAL "")
    message (FATAL_ERROR "In test ${TEST_NAME} NUM_CPUS not defined. Quitting")
  endif ()

  homme_add_system_test (${TEST_NAME}
    EXEC          ${EXEC_NAME}
    RANKS         ${NUM_CPUS}
    MAX_RANKS     ${MAX_NUM_CPUS}
    NAMELISTS     ${NAMELIST_FILES}
    VCOORD_FILES  ${VCOORD_FILES}
    MESH_FILES    ${MESH_FILES}
    NCL_FILES     ${NCL_FILES}
    OUTPUTS       ${NC_OUTPUT_FILES}
    OMP_NAMELISTS ${OMP_NAMELIST_FILES}
    OMP_THREADS   ${OMP_NUM_THREADS}
    REF_OUTPUTS   ${NC_OUTPUT_REF}
    CHECK_OUTPUTS ${NC_OUTPUT_CHECKREF}
    REF_TOL       ${TESTCASE_REF_TOL}
    LABELS        ${THIS_TEST_LABELS})

  # Error norms checks requested by the test file (via check_transport_error_norms)
  if (TARGET ${EXEC_NAME})
    foreach (args IN LISTS THIS_TEST_ERROR_NORM_CHECKS)
      string (REPLACE "|" ";" args "${args}")
      homme_add_error_norms_test (${args})
    endforeach()
  endif()
endmacro()

# Called by test description files: check the error norms printed by the test
macro (check_transport_error_norms TEST_NAME TCEN_ERROR_ANCHOR TCEN_FILENAME TCEN_UPPER_BOUNDS)
  # The test may not be created (if its exec is not built), so just record the request
  string (REPLACE ";" "," _bounds "${TCEN_UPPER_BOUNDS}")
  list (APPEND THIS_TEST_ERROR_NORM_CHECKS "${TEST_NAME}|${TCEN_ERROR_ANCHOR}|${TCEN_FILENAME}|${_bounds}")
endmacro()

macro (createTests testList)
  foreach (test ${${testList}})
    createTest(${test})
  endforeach ()
endmacro ()

macro (createTestsWithProfile testList testsProfile)
  set (PROFILE ${testsProfile})
  foreach (test ${${testList}})
    createTest(${test})
  endforeach ()
  unset (PROFILE)
endmacro ()

macro (resetTestVariables)
  foreach (var IN ITEMS TEST_NAME EXEC_NAME NAMELIST_FILES VCOORD_FILES MESH_FILES NCL_FILES
                        NC_OUTPUT_FILES NC_OUTPUT_REF NC_OUTPUT_CHECKREF TESTCASE_REF_TOL NUM_CPUS
                        MAX_NUM_CPUS OMP_NAMELIST_FILES OMP_NUM_THREADS THIS_TEST_ERROR_NORM_CHECKS)
    set (${var})
  endforeach()
endmacro ()

# Customize tests lengths/resolution, depending on the testing profile
macro (set_homme_tests_parameters testFile profile)
  # theta tests for kokkos targets
  # HV scaled at 3.0 rate:
  # ne30 1e15, ne2 3.4e18, ne4 4.2e17
  # ne6 1.25e17, ne13 1.23e16, ne1024 2.5e10
  if ("${testFile}" MATCHES "theta-f")
    if ("${profile}" STREQUAL "dev")
      set (HOMME_TEST_NE 2)
      set (HOMME_TEST_NDAYS 1)
      set (HOMME_TEST_NU 3.4e18)
    elseif ("${profile}" STREQUAL "short")
      set (HOMME_TEST_NE 6)
      set (HOMME_TEST_NDAYS 1) #should be 9 but for GB reduce
      set (HOMME_TEST_NU 1.25e17)
    else ()
      set (HOMME_TEST_NE 13)
      set (HOMME_TEST_NDAYS 1) #should be 9 but for GB reduce
      set (HOMME_TEST_NU 1.2e16)
    endif ()
  # preqx tests for kokkos targets
  elseif ("${testFile}" MATCHES "preqx-nlev")
    if ("${testFile}" MATCHES ".*-tensorhv-*")
      if ("${profile}" STREQUAL "dev")
        set (HOMME_TEST_NE 2)
        set (HOMME_TEST_NDAYS 1)
      elseif ("${profile}" STREQUAL "short")
        set (HOMME_TEST_NE 4)
        set (HOMME_TEST_NDAYS 3)
      else ()
        set (HOMME_TEST_NE 12)
        set (HOMME_TEST_NDAYS 3)
      endif ()
    else ()
      if ("${profile}" STREQUAL "dev")
        set (HOMME_TEST_NE 2)
        set (HOMME_TEST_NDAYS 1)
      elseif ("${profile}" STREQUAL "short")
        set (HOMME_TEST_NE 4)
        set (HOMME_TEST_NDAYS 9)
      else ()
        set (HOMME_TEST_NE 12)
        set (HOMME_TEST_NDAYS 9)
      endif ()
    endif()
  endif ()
endmacro ()

# Make a list of all testing profiles no more intensive than the given profile.
function (make_profiles_up_to profile profiles)
  string (TOLOWER "${profile}" profile_ci)
  if ("${profile_ci}" STREQUAL "dev")
    set (tmp dev)
  elseif ("${profile_ci}" STREQUAL "short")
    set (tmp dev short)
  elseif ("${profile_ci}" STREQUAL "nightly")
    set (tmp dev short nightly)
  else ()
    message (FATAL_ERROR "Testing profile '${profile}' not implemented.")
  endif ()
  set (${profiles} "${tmp}" PARENT_SCOPE)
endfunction ()

# For each test X in the list, create a cxx-vs-f90 comparison test between X and X-kokkos
macro (CREATE_CXX_VS_F90_TESTS_WITH_PROFILE TESTS_LIST testProfile)
  foreach (TEST ${${TESTS_LIST}})
    set_homme_tests_parameters(${TEST} ${testProfile})
    set (PROFILE ${testProfile})
    resetTestVariables()
    set (HOMME_ROOT ${HOMME_SOURCE_DIR})
    include (${HOMME_SOURCE_DIR}/test/reg_test/run_tests/${TEST}.cmake)

    set (TEST_NAME_SUFFIX)
    if ("${TEST}" MATCHES "theta-f")
      set (TEST_NAME_SUFFIX "-ne${HOMME_TEST_NE}-nu${HOMME_TEST_NU}-ndays${HOMME_TEST_NDAYS}")
    elseif ("${TEST}" MATCHES "preqx-nlev")
      set (TEST_NAME_SUFFIX "-ne${HOMME_TEST_NE}-ndays${HOMME_TEST_NDAYS}")
    endif ()

    message (STATUS "Creating cxx-f90 comparison test ${TEST}${TEST_NAME_SUFFIX}_cxx_vs_f90")
    homme_add_cxx_vs_f90_test ("${TEST}${TEST_NAME_SUFFIX}_cxx_vs_f90"
      F90_TEST "${TEST}${TEST_NAME_SUFFIX}"
      CXX_TEST "${TEST}-kokkos${TEST_NAME_SUFFIX}"
      OUTPUTS  ${NC_OUTPUT_FILES}
      LABELS   ${testProfile})
  endforeach ()
  unset (PROFILE)
endmacro ()

# Sets NUM_CPUS = min(USE_NUM_PROCS, MPI_MAX_PROCS), used by the unit tests.
# MPI_MAX_PROCS: maximum MPI ranks the test supports (default: 24).
macro (set_mpi_num_cpus)
  if (NOT DEFINED MPI_MAX_PROCS)
    set (MPI_MAX_PROCS 24)
  endif()
  if (USE_NUM_PROCS)
    if (${USE_NUM_PROCS} GREATER ${MPI_MAX_PROCS})
      message (STATUS "Test ${TEST_NAME}: lowering NUM_CPUS from ${USE_NUM_PROCS} to ${MPI_MAX_PROCS} (MPI_MAX_PROCS=${MPI_MAX_PROCS})")
      set (NUM_CPUS ${MPI_MAX_PROCS})
    else ()
      set (NUM_CPUS ${USE_NUM_PROCS})
    endif ()
  else ()
    set (NUM_CPUS ${MPI_MAX_PROCS})
  endif ()
endmacro()
