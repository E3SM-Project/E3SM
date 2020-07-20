set (MPICXX_WRAPPER_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/../bin)
macro (EkatSetNvccWrapper)
  SetMpiCxxBackendCompilerVarName("MPI_CXX_BACKEND_COMPILER_VAR_NAME")
  # Before starting the project, wrap mpicxx in the scream_mpicxx script, which
  # takes care of setting OMPI_CXX to point to nvcc_wrapper
  configure_file(${MPICXX_WRAPPER_SOURCE_DIR}/ekat_mpicxx.in ${CMAKE_BINARY_DIR}/bin/ekat_mpicxx @ONLY)
  set(CMAKE_CXX_COMPILER ${CMAKE_BINARY_DIR}/bin/ekat_mpicxx CACHE STRING "" FORCE)
endmacro()
