# This function takes the following arguments:
#    - target_name: the name of the executable
#    - target_srcs: a list of src files for the executable
#      Note: no need to include catch_main; this macro will add it
#    - scream_libs: a list of scream libraries needed by the executable
#    - config_defines: a list of additional defines for the compiler
#    - num_mpi_ranks: the number of mpi ranks for the test
#    - config_defines: specific define directives to be passed to the compiler

FUNCTION(CreateUnitTest target_name target_srcs scream_libs num_mpi_ranks config_defines)
  # Set link directories (must be done BEFORE add_executable is called)
  link_directories(${SCREAM_TPL_LIBRARY_DIRS})

  # Create the executable
  add_executable (${target_name} ${target_srcs} ${SCREAM_SRC_DIR}/share/util/scream_catch_main.cpp)

  # Set all target properties
  target_link_libraries(${target_name} ${scream_libs} ${SCREAM_TPL_LIBRARIES})
  target_include_directories(${target_name} PUBLIC ${SCREAM_INCLUDE_DIRS} ${CATCH_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
  set_target_properties(${target_name} PROPERTIES LINK_FLAGS "${SCREAM_LINK_FLAGS}")
  set_target_properties(${target_name} PROPERTIES Fortran_MODULE_DIRECTORY ${SCREAM_F90_MODULES})
  set_target_properties(${target_name} PROPERTIES COMPILE_DEFINITIONS "${config_defines}")

  # Create the test
  IF (${num_mpi_ranks} GREATER 1)
    ADD_TEST(${target_name}_ut mpiexec -np ${num_mpi_ranks} ${SCREAM_MPI_EXTRA_ARGS} ${target_name})
  ELSE()
    ADD_TEST(${target_name}_ut ${target_name})
  ENDIF()

ENDFUNCTION(CreateUnitTest)
