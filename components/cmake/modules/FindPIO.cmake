# - Try to find PIO
#
# Placeholder code until spio gets a proper config.cmake
#
# Once done, this will define:
#
#   The "spio" target
#

if (TARGET spio)
  return()
endif()

if (NOT PIO_LIBDIR)
  set(PIO_LIBDIR "${INSTALL_SHAREDPATH}/lib")
endif()

if (PIO_VERSION STREQUAL 2)
  # This is a pio2 library
  set(PIOLIBS "${PIO_LIBDIR}/libpiof.a;${PIO_LIBDIR}/libpioc.a")
else()
  # This is a pio1 library
  set(PIOLIBS "${PIO_LIBDIR}/libpio.a")
endif()

# Handle gptl. Just hardcode it for now.
list(APPEND PIOLIBS "${INSTALL_SHAREDPATH}/lib/libgptl.a")

find_package(NETCDF REQUIRED)

list(APPEND PIOLIBS netcdf)

if (MPILIB STREQUAL "mpi-serial")
  list(APPEND PIOLIBS "${INSTALL_SHAREDPATH}/lib/libmpi-serial.a")
else()
  find_package(MPI REQUIRED COMPONENTS C Fortran)
  list(APPEND PIOLIBS MPI::MPI_C MPI::MPI_Fortran)
endif()

# Create the interface library, and set target properties
add_library(spio INTERFACE)
target_link_libraries(spio INTERFACE ${PIOLIBS})
