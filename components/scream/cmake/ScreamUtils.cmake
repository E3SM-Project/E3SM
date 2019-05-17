include(CMakeParseArguments) # Needed for backwards compatibility

# This function takes the following arguments:
#    - target_name: the name of the executable
#    - target_srcs: a list of src files for the executable
#      Note: no need to include catch_main; this macro will add it
#    - scream_libs: a list of scream libraries needed by the executable
#    - config_defs (optional): a list of additional defines for the compiler, default is nothing
#    - mpi_ranks (optional): the number of mpi ranks for the test, if 2 values, it's a range, if 3, it's a range plus increment. default is np=1
#    - threads (optional): the number of threads for the test, if 2 values, it's a range, if 3, it's a range plus an increment. default is 1 thread
#      Note: One test will be created per combination of valid mpi-rank and thread value

FUNCTION(CreateUnitTest target_name target_srcs scream_libs)
  set(options OPTIONAL EXCLUDE_MAIN_CPP)
  set(oneValueArgs NAMELIST_FILE)
  set(multiValueArgs MPI_RANKS THREADS CONFIG_DEFS INCLUDE_DIRS)
  cmake_parse_arguments(CreateUnitTest "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # Set link directories (must be done BEFORE add_executable is called)
  link_directories(${SCREAM_TPL_LIBRARY_DIRS})

  # Create the executable
  set (SRC_MAIN ${SCREAM_SRC_DIR}/share/util/scream_catch_main.cpp)
  if (CreateUnitTest_EXCLUDE_MAIN_CPP)
    set (SRC_MAIN)
  endif ()
  add_executable (${target_name} ${target_srcs} ${SRC_MAIN})

  SET (TEST_INCLUDE_DIRS
       ${SCREAM_INCLUDE_DIRS}
       ${CATCH_INCLUDE_DIR}
       ${CMAKE_CURRENT_SOURCE_DIR}
       ${CMAKE_CURRENT_BINARY_DIR}
       ${SCREAM_F90_MODULES}
       ${CreateUnitTest_INCLUDE_DIRS}
  )

  # Set all target properties
  target_link_libraries(${target_name} ${scream_libs} ${SCREAM_TPL_LIBRARIES})
  target_include_directories(${target_name} PUBLIC ${TEST_INCLUDE_DIRS})
  set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${SCREAM_LINK_FLAGS}")
  set_target_properties(${target_name} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${target_name}_modules)
  if (CreateUnitTest_CONFIG_DEFS)
    set_target_properties(${target_name} PROPERTIES COMPILE_DEFINITIONS "${CreateUnitTest_CONFIG_DEFS}")
  endif()

  list(LENGTH CreateUnitTest_MPI_RANKS NUM_MPI_RANK_ARGS)
  list(LENGTH CreateUnitTest_THREADS   NUM_THREAD_ARGS)

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
    list(GET CreateUnitTest_MPI_RANKS 0 RETURN_VAL)
    set(MPI_START_RANK ${RETURN_VAL})
    set(MPI_END_RANK ${RETURN_VAL})
  elseif(NUM_MPI_RANK_ARGS EQUAL 2)
    list(GET CreateUnitTest_MPI_RANKS 0 RETURN_VAL)
    set(MPI_START_RANK ${RETURN_VAL})
    list(GET CreateUnitTest_MPI_RANKS 1 RETURN_VAL)
    set(MPI_END_RANK ${RETURN_VAL})
  else()
    list(GET CreateUnitTest_MPI_RANKS 0 RETURN_VAL)
    set(MPI_START_RANK ${RETURN_VAL})
    list(GET CreateUnitTest_MPI_RANKS 1 RETURN_VAL)
    set(MPI_END_RANK ${RETURN_VAL})
    list(GET CreateUnitTest_MPI_RANKS 2 RETURN_VAL)
    set(MPI_INCREMENT ${RETURN_VAL})
  endif()

  if (NUM_THREAD_ARGS EQUAL 0)
  elseif(NUM_THREAD_ARGS EQUAL 1)
    list(GET CreateUnitTest_THREADS 0 RETURN_VAL)
    set(THREAD_START ${RETURN_VAL})
    set(THREAD_END ${RETURN_VAL})
  elseif(NUM_THREAD_ARGS EQUAL 2)
    list(GET CreateUnitTest_THREADS 0 RETURN_VAL)
    set(THREAD_START ${RETURN_VAL})
    list(GET CreateUnitTest_THREADS 1 RETURN_VAL)
    set(THREAD_END ${RETURN_VAL})
  else()
    list(GET CreateUnitTest_THREADS 0 RETURN_VAL)
    set(THREAD_START ${RETURN_VAL})
    list(GET CreateUnitTest_THREADS 1 RETURN_VAL)
    set(THREAD_END ${RETURN_VAL})
    list(GET CreateUnitTest_THREADS 2 RETURN_VAL)
    set(THREAD_INCREMENT ${RETURN_VAL})
  endif()

  set(CURR_RANKS ${MPI_START_RANK})
  set(CURR_THREADS ${THREAD_START})

  IF (CreateUnitTest_NAMELIST_FILE)
    SET  (invokeExec "./${target_name} < ${CreateUnitTest_NAMELIST_FILE}")
  ELSE()
    SET  (invokeExec "./${target_name}")
  ENDIF()

  while (NOT CURR_RANKS GREATER ${MPI_END_RANK})
    while (NOT CURR_THREADS GREATER ${THREAD_END})
      # Create the test
      IF (${CURR_RANKS} GREATER 1)
        ADD_TEST(NAME ${target_name}_ut_np${CURR_RANKS}_omp${CURR_THREADS}
                 COMMAND sh -c "mpiexec -np ${CURR_RANKS} ${SCREAM_MPI_EXTRA_ARGS} ${invokeExec}")
      ELSE()
        ADD_TEST(NAME ${target_name}_ut_np1_omp${CURR_THREADS}
                 COMMAND sh -c "${invokeExec}")
      ENDIF()
      set_tests_properties(${target_name}_ut_np${CURR_RANKS}_omp${CURR_THREADS} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=${CURR_THREADS})
      MATH(EXPR CURR_THREADS "${CURR_THREADS}+${THREAD_INCREMENT}")
    endwhile()
    set(CURR_THREADS ${THREAD_START})
    MATH(EXPR CURR_RANKS "${CURR_RANKS}+${MPI_INCREMENT}")
  endwhile()

ENDFUNCTION(CreateUnitTest)
