

find_path(Netcdf_INCLUDE_DIR 
          netcdf.h
          #PATHS ${Homme_NETCDF_DIR} ${NETCDF_DIR}
          PATHS ${NETCDF_DIR}
          PATH_SUFFIXES include
          NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

find_library(Netcdf_LIBRARY 
             NAMES libnetcdf.a netcdf
             HINTS ${Netcdf_INCLUDE_DIR}/../lib)

find_library(NetcdfF_LIBRARY 
             NAMES libnetcdff.a netcdff
             HINTS ${Netcdf_INCLUDE_DIR}/../lib)

set(Netcdf_LIBRARIES ${NetcdfF_LIBRARY} ${Netcdf_LIBRARY})

# Now check to see if Netcdf can be linked to without any other (static) libraries
# This check creates a small example program and tries to link with ncopen
#   it handles some symbol mangling
INCLUDE(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(${Netcdf_LIBRARY} ncopen "" NETCDF_SELF_CONTAINED)

IF ("${NETCDF_SELF_CONTAINED}")
  SET(NETCDF_REQUIRE_HDF5 FALSE)
ELSE ()
  MESSAGE(STATUS "This version of Netcdf requires HDF5")
  SET(NETCDF_REQUIRE_HDF5 TRUE)
ENDIF ()

IF (${NETCDF_REQUIRE_HDF5}) 

  find_path(HDF5_INCLUDE_DIR
            hdf5.h
            PATHS ${Homme_HDF5_DIR} ${HDF5_DIR}
            HINTS ${Homme_Hint_Paths}
            PATH_SUFFIXES include
            NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  find_library(HDF5_LIBRARY
               NAMES libhdf5.a hdf5
               PATHS ${HDF5_DIR}
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

  find_library(HDF5hl_LIBRARY
               NAMES libhdf5_hl.a hdf5_hl
               PATHS ${HDF5_DIR}
               PATH_SUFFIXES lib lib64
               NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)


  if(${HDF5_LIBRARY} STREQUAL "HDF5_LIBRARY-NOTFOUND" OR ${HDF5hl_LIBRARY} STREQUAL "HDF5hl_LIBRARY-NOTFOUND")
    set(HDF5_FOUND OFF)
    MESSAGE(FATAL_ERROR "HDF5 not found")
  else()
    set(HDF5_FOUND ON)
    MESSAGE(STATUS "Found HDF5: ${HDF5hl_LIBRARY} ${HDF5_LIBRARY}")
    set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${HDF5hl_LIBRARY} ${HDF5_LIBRARY})
  endif()

  # Check to see which dependencies (ZLIB, SZIP) hdf5 has
  INCLUDE(CheckSymbolExists)
  CHECK_SYMBOL_EXISTS(H5_HAVE_SZLIB_H "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_REQUIRE_SZ)
  CHECK_SYMBOL_EXISTS(H5_HAVE_ZLIB_H "${HDF5_INCLUDE_DIR}/H5pubconf.h" HDF5_REQUIRE_ZLIB)

  #MESSAGE(STATUS "HDF5_REQUIRE_SZ=${HDF5_REQUIRE_SZ}")
  #MESSAGE(STATUS "HDF5_REQUIRE_ZLIB=${HDF5_REQUIRE_ZLIB}")

  IF (${HDF5_REQUIRE_ZLIB})

    MESSAGE(STATUS "HDF5 Requires ZLIB")

    find_library(ZLIB_LIBRARY
                 NAMES libz.a z
                 PATHS ${ZLIB_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
    if(${ZLIB_LIBRARY} STREQUAL "ZLIB_LIBRARY-NOTFOUND")
      set(ZLIB_FOUND OFF)
      MESSAGE(FATAL_ERROR "ZLIB Not found")
    else()
      set(ZLIB_FOUND ON)
      MESSAGE(STATUS "Found ZLIB: ${ZLIB_LIBRARY}")
      set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${ZLIB_LIBRARY})
    endif()

  ENDIF ()

  IF (${HDF5_REQUIRE_SZ})

    MESSAGE(STATUS "HDF5 Requires SZIP")

    find_library(SZIP_LIBRARY
                 NAMES libsz.a sz
                 PATHS ${SZIP_DIR}
                 PATH_SUFFIXES lib lib64
                 NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)
    if(${SZIP_LIBRARY} STREQUAL "SZIP_LIBRARY-NOTFOUND")
      set(ZLIB_FOUND OFF)
      MESSAGE(FATAL_ERROR "SZIP Not found")
    else()
      set(ZLIB_FOUND ON)
      MESSAGE(STATUS "Found SZIP: ${SZIP_LIBRARY}")
      set(Netcdf_LIBRARIES ${Netcdf_LIBRARIES} ${SZIP_LIBRARY})
    endif()

  ENDIF ()

  # Need to handle include dirs
  set(Netcdf_INCLUDE_DIRS ${Netcdf_INCLUDE_DIR})

ENDIF()

set(Netcdf_FOUND ON)
MESSAGE(STATUS "Found Netcdf: ${Netcdf_LIBRARIES}")

IF(Netcdf_FIND_REQUIRED AND NOT Netcdf_FOUND)
  MESSAGE(FATAL_ERROR "Did not find required library Netcdf")
ELSEIF(Netcdf_FIND_REQUIRED AND NOT HDF5_FOUND)
  MESSAGE(STATUS "Found Netcdf and all dependencies")
ENDIF()
