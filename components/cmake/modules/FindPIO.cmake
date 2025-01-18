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
  if (DEFINED ENV{ADIOS2_ROOT})
    list(APPEND PIOLIBS "${PIO_LIBDIR}/libadios2pio-nm-lib.a")
  endif()
else()
  # This is a pio1 library
  set(PIOLIBS "${PIO_LIBDIR}/libpio.a")
endif()

# Handle gptl. Just hardcode it for now.
list(APPEND PIOLIBS "${INSTALL_SHAREDPATH}/lib/libgptl.a")

find_package(NETCDF REQUIRED)
# Check if scorpio has hdf5 enabled
if (DEFINED ENV{HDF5_ROOT})
  if (DEFINED ENV{HDF5_USE_STATIC_LIBRARIES})
    set(HDF5_USE_STATIC_LIBRARIES On)
  endif()
  find_package(HDF5 REQUIRED COMPONENTS C HL)
endif()

# Not all machines/PIO installations use ADIOS but, for now,
# we can assume that an MPI case with ADIOS2_ROOT set is probably
# using adios.
if (NOT MPILIB STREQUAL "mpi-serial" AND DEFINED ENV{ADIOS2_ROOT})
  if(DEFINED ENV{BLOSC2_ROOT})
    set(ENV{Blosc2_DIR} "$ENV{BLOSC2_ROOT}")
  endif()
  find_package(MPI REQUIRED COMPONENTS C)
  find_package(ADIOS2 REQUIRED COMPONENTS C)
  list(APPEND PIOLIBS adios2::adios2)
endif()

list(APPEND PIOLIBS netcdf)

if (MPILIB STREQUAL "mpi-serial")
  list(APPEND PIOLIBS "${INSTALL_SHAREDPATH}/lib/libmpi-serial.a")
else()
  find_package(MPI REQUIRED COMPONENTS C Fortran)
  list(APPEND PIOLIBS MPI::MPI_C MPI::MPI_Fortran)
endif()

if (DEFINED ENV{HDF5_ROOT})
  list(APPEND PIOLIBS ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
endif()

# Create the interface library, and set target properties
add_library(spio INTERFACE)
target_link_libraries(spio INTERFACE ${PIOLIBS})
