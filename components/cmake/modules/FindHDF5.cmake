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

if (NOT EXISTS "${HDF5_ROOT}/lib" AND NOT EXISTS "${HDF5_ROOT}/lib64")
  message(FATAL_ERROR "HDF5_ROOT does not contain a lib or lib64 directory")
endif()

# Preserve original CMake find library suffixes
set(ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES "${CMAKE_FIND_LIBRARY_SUFFIXES}")

if (HDF5_USE_STATIC_LIBRARIES)
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
else()
  set(CMAKE_FIND_LIBRARY_SUFFIXES .so .a)
endif()

find_library(HDF5_LIBRARIES NAMES hdf5 HINTS "${HDF5_ROOT}/lib" "${HDF5_ROOT}/lib64")
find_library(HDF5_HL_LIBRARIES NAMES hdf5_hl HINTS "${HDF5_ROOT}/lib" "${HDF5_ROOT}/lib64")

# Restore original CMake find library suffixes
set(CMAKE_FIND_LIBRARY_SUFFIXES "${ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES}")

if (NOT EXISTS "${HDF5_ROOT}/include")
  message(FATAL_ERROR "HDF5_ROOT does not contain an include directory")
endif()

find_path(HDF5_INCLUDE_DIR hdf5.h HINTS "${HDF5_ROOT}/include")

# Create the interface library, and set target properties
# For static libraries, link with HDF5_HL_LIBRARIES before HDF5_LIBRARIES
add_library(hdf5 INTERFACE)
target_link_libraries(hdf5 INTERFACE ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
target_include_directories(hdf5 INTERFACE "${HDF5_INCLUDE_DIR}")
