# - Try to find HDF5
#
# CMake's built-in FindHDF5 module may fail with certain compilers, so
# we provide a custom module that avoids relying on compiler wrappers.
#
# Once done, this will define:
#
#   The "hdf5" target
#

if (TARGET hdf5)
  return()
endif()

set(HDF5_ROOT $ENV{HDF5_ROOT})

# Verify HDF5 installation by checking for h5dump
if (NOT EXISTS "${HDF5_ROOT}/bin/h5dump")
  message(FATAL_ERROR "h5dump not found in ${HDF5_ROOT}/bin")
endif()

# Include and library paths
find_path(HDF5_INCLUDE_DIR hdf5.h HINTS ${HDF5_ROOT}/include)
find_library(HDF5_LIBRARIES NAMES hdf5 HINTS ${HDF5_ROOT}/lib ${HDF5_ROOT}/lib64)
find_library(HDF5_HL_LIBRARIES NAMES hdf5_hl HINTS ${HDF5_ROOT}/lib ${HDF5_ROOT}/lib64)

# Create the interface library, and set target properties
add_library(hdf5 INTERFACE)
target_include_directories(hdf5 INTERFACE ${HDF5_INCLUDE_DIR})
target_link_libraries(hdf5 INTERFACE ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
