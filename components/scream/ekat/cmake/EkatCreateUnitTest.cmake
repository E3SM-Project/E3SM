include(CMakeParseArguments) # Needed for backwards compatibility

# This function takes the following mandatory arguments:
#    - target_name: the name of the executable
#    - target_srcs: a list of src files for the executable
#      Note: no need to include ekat_catch_main; this macro will add it
# and the following optional arguments (to be passed as ARG_NAME "ARG_VAL")
#    - MPI_RANKS: the number of mpi ranks for the test.
#      Note: if 2 values, it's a range, if 3, it's a range plus increment. default is np=1
#    - THREADS: the number of threads for the test
#      Note: if 2 values, it's a range, if 3, it's a range plus an increment. default is 1 thread
#      Note: for each combination of valid mpi-rank and thread value, a new test will be created,
#            with suffix '_npN_omp_M', with N numver of mpi ranks, and M number of omp threads.
#    - MPI_EXEC_NAME: name of the mpi launcher (usually, mpiexe or mpirun, but may be another wrapper)
#    - MPI_NP_FLAG: the flag used to specify the number of mpi ranks (usually, -np or -n)
#    - MPI_EXTRA_ARGS: additional args to be forwarded to the mpi launches (e.g., --map-by, --bind-to, ...)
#    - COMPILE_DEFS: a list of additional defines for the compiler
#    - COMPILER_FLAGS: a list of additional flags for the compiler
#    - LIBS: a list of libraries needed by the executable
#    - LIBS_DIRS: a list of directories to add to the linker search path
#    - LINKER_FLAGS: a list of additional flags for the linker
#    - LABELS: a set of labels to attach to the test

# Note: we hace to set this variable here, so CMAKE_CURRENT_LIST_DIR gets the
#       directory of this file. If we did it inside the function, it would get
#       the directory from where the function is called
set (CATCH_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/../extern/catch2/include)

function(EkatCreateUnitTest target_name target_srcs)
  set(options OPTIONAL EXCLUDE_MAIN_CPP)
  set(oneValueArgs DEP MPI_EXEC_NAME MPI_NP_FLAG)
  set(multiValueArgs
    MPI_RANKS THREADS
    MPI_EXTRA_ARGS EXE_ARGS 
    COMPILER_DEFS INCLUDE_DIRS COMPILER_FLAGS
    LIBS LIBS_DIRS LINKER_FLAGS
    LABELS)

  # ecut = Ekat Create Unit Test
  cmake_parse_arguments(ecut "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # Set link directories (must be done BEFORE add_executable is called)
  if (ecut_LIBS_DIRS)
    link_directories(${ecut_LIBS_DIRS})
  endif()

  # Create the executable
  if (ecut_EXCLUDE_MAIN_CPP)
    set (SRC_MAIN)
  else()
    set (SRC_MAIN ${EKAT_SOURCE_DIR}/src/ekat/util/ekat_catch_main.cpp)
  endif ()
  add_executable (${target_name} ${target_srcs} ${SRC_MAIN})

  set (TEST_INCLUDE_DIRS
       ${CATCH_INCLUDE_DIR}
       ${CMAKE_CURRENT_SOURCE_DIR}
       ${CMAKE_CURRENT_BINARY_DIR}
       ${ecut_INCLUDE_DIRS}
  )

  # Set all target properties
  target_include_directories(${target_name} PUBLIC ${TEST_INCLUDE_DIRS})
  if (ecut_LIBS)
    target_link_libraries(${target_name} "${ecut_LIBS}")
  endif()
  if (ecut_LINKER_FLAGS)
    set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${ecut_LINKER_FLAGS}")
  endif()
  set_target_properties(${target_name} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${target_name}_modules)
  if (ecut_COMPILER_DEFS)
    set_target_properties(${target_name} PROPERTIES COMPILE_DEFINITIONS "${ecut_COMPILER_DEFS}")
  endif()

  list(LENGTH ecut_MPI_RANKS NUM_MPI_RANK_ARGS)
  list(LENGTH ecut_THREADS   NUM_THREAD_ARGS)

  if (NUM_MPI_RANK_ARGS GREATER 3)
    message(FATAL_ERROR "Too many mpi arguments for ${target_name}")
  endif()
  if (NUM_THREAD_ARGS GREATER 3)
    message(FATAL_ERROR "Too many thread arguments for ${target_name}")
  endif()

  set(MPI_START_RANK 1)
  set(MPI_END_RANK 1)
  set(MPI_INCREMENT 1)

  set(THREAD_START 1)
  set(THREAD_END 1)
  set(THREAD_INCREMENT 1)

  if (NUM_MPI_RANK_ARGS EQUAL 0)
  elseif(NUM_MPI_RANK_ARGS EQUAL 1)
    list(GET ecut_MPI_RANKS 0 RETURN_VAL)
    set(MPI_START_RANK ${RETURN_VAL})
    set(MPI_END_RANK ${RETURN_VAL})
  elseif(NUM_MPI_RANK_ARGS EQUAL 2)
    list(GET ecut_MPI_RANKS 0 RETURN_VAL)
    set(MPI_START_RANK ${RETURN_VAL})
    list(GET ecut_MPI_RANKS 1 RETURN_VAL)
    set(MPI_END_RANK ${RETURN_VAL})
  else()
    list(GET ecut_MPI_RANKS 0 RETURN_VAL)
    set(MPI_START_RANK ${RETURN_VAL})
    list(GET ecut_MPI_RANKS 1 RETURN_VAL)
    set(MPI_END_RANK ${RETURN_VAL})
    list(GET ecut_MPI_RANKS 2 RETURN_VAL)
    set(MPI_INCREMENT ${RETURN_VAL})
  endif()

  if (NUM_THREAD_ARGS EQUAL 0)
  elseif(NUM_THREAD_ARGS EQUAL 1)
    list(GET ecut_THREADS 0 RETURN_VAL)
    set(THREAD_START ${RETURN_VAL})
    set(THREAD_END ${RETURN_VAL})
  elseif(NUM_THREAD_ARGS EQUAL 2)
    list(GET ecut_THREADS 0 RETURN_VAL)
    set(THREAD_START ${RETURN_VAL})
    list(GET ecut_THREADS 1 RETURN_VAL)
    set(THREAD_END ${RETURN_VAL})
  else()
    list(GET ecut_THREADS 0 RETURN_VAL)
    set(THREAD_START ${RETURN_VAL})
    list(GET ecut_THREADS 1 RETURN_VAL)
    set(THREAD_END ${RETURN_VAL})
    list(GET ecut_THREADS 2 RETURN_VAL)
    set(THREAD_INCREMENT ${RETURN_VAL})
  endif()

  if (NOT MPI_START_RANK GREATER 0)
    message (FATAL_ERROR "Error! MPI_START_RANK is <=0.")
  endif()
  if (NOT MPI_END_RANK GREATER 0)
    message (FATAL_ERROR "Error! MPI_END_RANK is <=0.")
  endif()
  if (MPI_INCREMENT GREATER 0 AND MPI_START_RANK GREATER MPI_END_RANK)
    message (FATAL_ERROR "Error! MPI_START_RANK > MPI_END_RANK, but the increment is positive.")
  endif()
  if (MPI_INCREMENT LESS 0 AND MPI_START_RANK LESS MPI_END_RANK)
    message (FATAL_ERROR "Error! MPI_START_RANK < MPI_END_RANK, but the increment is negative.")
  endif()
  if (NOT THREAD_START GREATER 0)
    message (FATAL_ERROR "Error! THREAD_START is <=0.")
  endif()
  if (NOT THREAD_END GREATER 0)
    message (FATAL_ERROR "Error! THREAD_END is <=0.")
  endif()
  if (THREAD_INCREMENT GREATER 0 AND THREAD_START GREATER THREAD_END)
    message (FATAL_ERROR "Error! THREAD_START > THREAD_END, but the increment is positive.")
  endif()
  if (THREAD_INCREMENT LESS 0 AND THREAD_START LESS THREAD_END)
    message (FATAL_ERROR "Error! THREAD_START < THREAD_END, but the increment is negative.")
  endif()

  # Check both, in case user has negative increment
  if (MPI_END_RANK GREATER 1 OR MPI_START_RANK GREATER 1)
    if ("${ecut_MPI_EXEC_NAME}" STREQUAL "")
      set (ecut_MPI_EXEC_NAME "mpiexec")
    endif()
    if ("${ecut_MPI_NP_FLAG}" STREQUAL "")
      set (ecut_MPI_NP_FLAG "-n")
    endif()
  endif()

  set(CURR_RANKS ${MPI_START_RANK})
  set(CURR_THREADS ${THREAD_START})

  if (ecut_EXE_ARGS)
    set(invokeExec "./${target_name} ${ecut_EXE_ARGS}")
  else()
    set(invokeExec "./${target_name}")
  endif()

  while (NOT CURR_RANKS GREATER ${MPI_END_RANK})
    while (NOT CURR_THREADS GREATER ${THREAD_END})
      # Create the test
      set(FULL_TEST_NAME ${target_name}_ut_np${CURR_RANKS}_omp${CURR_THREADS})
      if (${CURR_RANKS} GREATER 1)
        add_test(NAME ${FULL_TEST_NAME}
                 COMMAND sh -c "${ecut_MPI_EXEC_NAME} ${ecut_MPI_NP_FLAG} ${CURR_RANKS} ${ecut_MPI_EXTRA_ARGS} ${invokeExec}")
      else()
        add_test(NAME ${FULL_TEST_NAME}
                 COMMAND sh -c "${invokeExec}")
      endif()
      math(EXPR CURR_CORES "${CURR_RANKS}*${CURR_THREADS}")
      set_tests_properties(${FULL_TEST_NAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=${CURR_THREADS} PROCESSORS ${CURR_CORES} PROCESSOR_AFFINITY True)
      if (ecut_DEP AND NOT ecut_DEP STREQUAL "${FULL_TEST_NAME}")
        set_tests_properties(${FULL_TEST_NAME} PROPERTIES DEPENDS ${ecut_DEP})
      endif()

      if (ecut_LABELS)
        set_tests_properties(${FULL_TEST_NAME} PROPERTIES LABELS "${ecut_LABELS}")
      endif()

      math(EXPR CURR_THREADS "${CURR_THREADS}+${THREAD_INCREMENT}")
    endwhile()
    set(CURR_THREADS ${THREAD_START})
    math(EXPR CURR_RANKS "${CURR_RANKS}+${MPI_INCREMENT}")
  endwhile()

endfunction(EkatCreateUnitTest)
