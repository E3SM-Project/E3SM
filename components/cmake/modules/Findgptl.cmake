# - Try to find gptl
#
# If gptl ever implements GptlConfig.cmake and installs as a proper
# cmake package, we will be able to remove this file

if (TARGET gptl)
  return()
endif()

# Find the build gptl library
if (NOT GPTL_LIBDIR)
  set(GPTL_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()
find_library(GPTL_LIB NAME gptl HINTS ${GPTL_LIBDIR})

# Create the interface target, and link the library
add_library (gptl INTERFACE)
target_link_libraries (gptl INTERFACE ${GTPL_LIB})
target_include_directories (gptl INTERFACE ${INSTALL_SHAREDPATH}/include) # for .mod files

# Link mpi or mpi-serial libs
if (MPILIB STREQUAL "mpi-serial")
  find_library(MPISERIAL_LIB NAME mpi-serial HINTS ${INSTALL_SHAREDPATH}/lib)
  target_link_libraries (gptl INTERFACE ${MPISERIAL_LIB})
else()
  find_package(MPI REQUIRED COMPONENTS C Fortran)
  target_link_libraries (gptl INTERFACE MPI::MPI_C MPI::MPI_Fortran)
endif()
