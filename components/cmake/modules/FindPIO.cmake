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

# Locate a library in PIO_LIBDIR (or the shared install dir), whether it was built static or shared
function(find_pio_library var name)
  find_library(${var} NAMES ${name} HINTS "${PIO_LIBDIR}" "${INSTALL_SHAREDPATH}/lib" NO_DEFAULT_PATH)
  if (NOT ${var})
    message(FATAL_ERROR "Could not find lib${name} in ${PIO_LIBDIR} or ${INSTALL_SHAREDPATH}/lib")
  endif()
endfunction()

if (PIO_VERSION STREQUAL 2)
  # This is a pio2 library
  find_pio_library(PIO_PIOF_LIBRARY piof)
  find_pio_library(PIO_PIOC_LIBRARY pioc)
  set(PIOLIBS ${PIO_PIOF_LIBRARY} ${PIO_PIOC_LIBRARY})
  if (DEFINED ENV{ADIOS2_ROOT})
    find_pio_library(PIO_ADIOS2PIO_LIBRARY adios2pio-nm-lib)
    list(APPEND PIOLIBS ${PIO_ADIOS2PIO_LIBRARY})
  endif()
else()
  # This is a pio1 library
  find_pio_library(PIO_PIO_LIBRARY pio)
  set(PIOLIBS ${PIO_PIO_LIBRARY})
endif()

# Handle gptl
find_pio_library(PIO_GPTL_LIBRARY gptl)
list(APPEND PIOLIBS ${PIO_GPTL_LIBRARY})

# Find the MPI library
if (NOT MPILIB STREQUAL "mpi-serial")
  find_package(MPI REQUIRED COMPONENTS C)
endif()

# ADIOS & HDF5 use Blosc2, ZFP for data compression
if (NOT MPILIB STREQUAL "mpi-serial" AND
      (DEFINED ENV{ADIOS2_ROOT} OR DEFINED ENV{HDF5_ROOT}))
  if(DEFINED ENV{BLOSC2_ROOT})
    set(ENV{Blosc2_DIR} "$ENV{BLOSC2_ROOT}")
  endif()
  if(DEFINED ENV{ZFP_ROOT})
    set(ENV{ZFP_DIR} "$ENV{ZFP_ROOT}")
  endif()
  if (DEFINED ENV{HDF5_ROOT)
    # Set paths to HDF5 filter libs for Blosc2, ZFP
    if(DEFINED ENV{H5Z_BLOSC2_ROOT})
      set(ENV{H5Z_BLOSC2_DIR} "$ENV{H5Z_BLOSC2_ROOT}")
    endif()
    if(DEFINED ENV{H5Z_ZFP_ROOT})
      set(ENV{H5Z_ZFP_DIR} "$ENV{H5Z_ZFP_ROOT}")
    endif()
  endif()
endif()

find_package(NETCDF REQUIRED)
# Check if scorpio has hdf5 enabled
if (DEFINED ENV{HDF5_ROOT})
  if (DEFINED ENV{HDF5_USE_STATIC_LIBRARIES})
    set(HDF5_USE_STATIC_LIBRARIES On)
  endif()
  find_package(HDF5 REQUIRED)

  # Find HDF5 filter libraries
  set(H5Z_ZFP_USE_STATIC_LIBS ON)
  list(APPEND CMAKE_PREFIX_PATH "$ENV{H5Z_ZFP_DIR}")
  find_package(H5Z_ZFP)

  set(H5Z_BLOSC2_USE_STATIC_LIBS ON)
  list(APPEND CMAKE_PREFIX_PATH "$ENV{H5Z_BLOSC2_DIR}")
  find_package(H5Z_BLOSC2)
endif()

# Not all machines/PIO installations use ADIOS but, for now,
# we can assume that an MPI case with ADIOS2_ROOT set is probably
# using adios.
if (NOT MPILIB STREQUAL "mpi-serial" AND DEFINED ENV{ADIOS2_ROOT})
  find_package(ADIOS2 REQUIRED COMPONENTS C)
  list(APPEND PIOLIBS adios2::adios2)
endif()

list(APPEND PIOLIBS netcdf)

if (MPILIB STREQUAL "mpi-serial")
  find_pio_library(PIO_MPI_SERIAL_LIBRARY mpi-serial)
  list(APPEND PIOLIBS ${PIO_MPI_SERIAL_LIBRARY})
else()
  find_package(MPI REQUIRED COMPONENTS C Fortran)
  list(APPEND PIOLIBS MPI::MPI_C MPI::MPI_Fortran)
endif()

if (DEFINED ENV{HDF5_ROOT})
  list(APPEND PIOLIBS hdf5)
  if (H5Z_ZFP_FOUND)
    list(APPEND PIOLIBS h5z_zfp::h5z_zfp)
  endif ()
  if (H5Z_BLOSC2_FOUND)
    list(APPEND PIOLIBS h5z_blosc2::h5z_blosc2)
  endif ()
endif()

# Create the interface library, and set target properties
add_library(spio INTERFACE)
target_link_libraries(spio INTERFACE ${PIOLIBS})
