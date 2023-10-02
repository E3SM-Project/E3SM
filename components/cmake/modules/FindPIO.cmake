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

find_package(NETCDF REQUIRED)

# Not all machines/PIO installations use ADIOS but, for now,
# we can assume that an MPI case with ADIOS2_ROOT set is probably
# using adios.
if (NOT MPILIB STREQUAL mpi-serial AND DEFINED $ENV{ADIOS2_ROOT}))
  find_package(MPI REQUIRED COMPONENTS C)
  find_package(ADIOS2 REQUIRED COMPONENTS C)
  list(APPEND PIOLIBS adios2::adios2)
endif()

list(APPEND PIOLIBS netcdf)

# Create the interface library, and set target properties
add_library(spio INTERFACE)
target_link_libraries(spio INTERFACE ${PIOLIBS})
