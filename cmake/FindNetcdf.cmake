MESSAGE(STATUS "Homme_NETCDF_DIR = ${Homme_NETCDF_DIR}")
find_path(Netcdf_INCLUDE_DIR 
          netcdf.h
          PATHS ${Homme_NETCDF_DIR} ${NETCDF_DIR}
          HINTS ${Homme_Hint_Paths}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

find_library(Netcdf_LIBRARY 
             NAMES libnetcdf.a netcdf
             HINTS ${Netcdf_INCLUDE_DIR}/../lib)

find_library(NetcdfF_LIBRARY 
             NAMES libnetcdff.a netcdff
             HINTS ${Netcdf_INCLUDE_DIR}/../lib)

find_library(HDF5_LIBRARY
             NAMES libhdf5.a hdf5
             PATHS ENV HDF5_DIR
             PATH_SUFFIXES lib lib64
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
find_library(HDF5hl_LIBRARY
             NAMES libhdf5_hl.a hdf5_hl
             PATHS ENV HDF5_DIR
             PATH_SUFFIXES lib lib64
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

if(${HDF5_LIBRARY} STREQUAL "HDF5_LIBRARY-NOTFOUND" OR ${HDF5hl_LIBRARY} STREQUAL "HDF5hl_LIBRARY-NOTFOUND")
  set(HDF5_FOUND OFF)
else()
  set(HDF5_FOUND ON)
  MESSAGE(STATUS "Found HDF5: ${HDF5hl_LIBRARY} ${HDF5_LIBRARY}")
endif()

find_library(ZLIB_LIBRARY
             NAMES libz.a z
             PATHS ENV Z_DIR
             PATH_SUFFIXES lib lib64
             NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
if(${ZLIB_LIBRARY} STREQUAL "ZLIB_LIBRARY-NOTFOUND")
  set(ZLIB_FOUND OFF)
else()
  set(ZLIB_FOUND ON)
  MESSAGE(STATUS "Found ZLIB: ${ZLIB_LIBRARY}")
endif()

if (${Netcdf_INCLUDE_DIR} STREQUAL "Netcdf_INCLUDE_DIR-NOTFOUND")
  set(Netcdf_FOUND OFF)
else()
  set(Netcdf_LIBRARIES ${NetcdfF_LIBRARY} ${Netcdf_LIBRARY})
  if (HDF5_FOUND)
    # Add in the hdf5 libraries
    set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${HDF5hl_LIBRARY} ${HDF5_LIBRARY})
  endif()
  if (ZLIB_FOUND)
    # Add in the zlib library
    set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${ZLIB_LIBRARY})
  endif() 
  set(Netcdf_INCLUDE_DIRS ${Netcdf_INCLUDE_DIR})
  set(Netcdf_FOUND ON)
  MESSAGE(STATUS "Found Netcdf: ${Netcdf_LIBRARIES}")
endif()

IF(Netcdf_FIND_REQUIRED AND NOT Netcdf_FOUND)
  MESSAGE(FATAL_ERROR "Did not find required library Netcdf")
ELSEIF(Netcdf_FIND_REQUIRED AND NOT HDF5_FOUND)
  MESSAGE(STATUS "Found Netcdf; did not find HDF5.")
ENDIF()
