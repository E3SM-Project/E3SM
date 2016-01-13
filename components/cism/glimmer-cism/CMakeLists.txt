# CMAKE File for CISM building against an installed Trilinos

cmake_minimum_required(VERSION 2.8.4)

OPTION (CISM_BUILD_CISM_DRIVER "Toggle to build cism_driver, on by default" ON)

OPTION (CISM_USE_TRILINOS "Toggle to use Trilinos: defaults to OFF" OFF)
OPTION (CISM_MPI_MODE "Toggle to Configure with MPI: defaults to ON" ON)
OPTION (CISM_SERIAL_MODE "Toggle to Configure in Serial mode: defaults to OFF " OFF)

OPTION (CISM_USE_MPI_WITH_SLAP "Toggle to use mpi when using SLAP solver, only relevant if CISM_SERIAL_MODE=ON: defaults to OFF" OFF)
OPTION (CISM_BUILD_SIMPLE_GLIDE "Toggle to build simple_glide, OFF by default" OFF)
OPTION (CISM_ENABLE_BISICLES "Toggle to build a BISICLES-capable cism_driver, off by default" OFF)

OPTION (CISM_BUILD_EXTRA_EXECUTABLES "Toggle to other executables, off by default" OFF)
OPTION (CISM_USE_GPTL_INSTRUMENTATION "Toggle to use GPTL instrumentation, on by default " ON)
OPTION (CISM_COUPLED "Toggle to build CISM for use with CESM, off by default" OFF)
OPTION (CISM_USE_DEFAULT_IO "Toggle to use default i/o files rather than running python script, off by default" OFF)
# OPTION (CISM_USE_CISM_FRONT_END "Toggle to use cism_driver or cism_cesm_interface with cism_front_end, off by default" OFF)

OPTION (CISM_STATIC_LINKING "Toggle to set static linking for executables, off by default" OFF)
OPTION (CISM_FORCE_FORTRAN_LINKER "Toggle to force using a fortran linker for building executables, off by default" OFF)
OPTION (CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES "Toggle to explicitly include the CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES on the link line, on by default" ON)

MESSAGE("CISM_USE_TRILINOS = ${CISM_USE_TRILINOS}")

# override CISM_SERIAL_MODE setting, if it conflicts with CISM_MPI_MODE setting:
IF (CISM_MPI_MODE)
  SET(CISM_SERIAL_MODE "OFF")
  SET(CISM_BUILD_EXTRA_EXECUTABLES "OFF")
ENDIF()

IF (NOT CISM_MPI_MODE AND NOT CISM_SERIAL_MODE)
  SET(CISM_MPI_MODE "ON")
ENDIF()

IF (CISM_BUILD_CISM_DRIVER)
  SET(CISM_USE_CISM_FRONT_END "ON")
ENDIF()

# set (or override) other options, if CISM_COUPLED is ON:
IF (CISM_COUPLED)
  SET(CISM_NO_EXECUTABLE "ON")
  SET(CISM_BUILD_SIMPLE_GLIDE "OFF")
  SET(CISM_ENABLE_BISICLES "OFF")
  SET(CISM_USE_CISM_FRONT_END "OFF")
  SET(CISM_USE_DEFAULT_IO "ON")
  ADD_DEFINITIONS(-DCCSMCOUPLED)
ENDIF()


IF (NOT DEFINED CISM_BINARY_DIR)
   SET(CISM_BINARY_DIR ${CMAKE_BINARY_DIR})
ENDIF()

MESSAGE("Building in: ${CISM_BINARY_DIR}")

IF (CISM_USE_TRILINOS)
  OPTION(CISM_TRILINOS_DIR "Required path to installed Trilinos")
  OPTION(CISM_NETCDF_DIR "Required path to installed Netcdf")


  IF (CISM_USE_GPTL_INSTRUMENTATION)
     IF (DEFINED CISM_TRILINOS_GPTL_DIR)
       SET(CISM_TRILINOS_DIR ${CISM_TRILINOS_GPTL_DIR})
     ENDIF()
  ENDIF()

  
  # Error check up front
  IF (NOT DEFINED CISM_TRILINOS_DIR)
    MESSAGE(FATAL_ERROR "\nCISM Error: cmake must define CISM_TRILINOS_DIR:
    (-D CISM_TRILINOS_DIR=<trilinos_install_prefix>)!")
  ENDIF()

  # Get Trilinos as one entity
  SET(CMAKE_PREFIX_PATH ${CISM_TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
  FIND_PACKAGE(Trilinos REQUIRED)

  IF (${Trilinos_VERSION} VERSION_LESS  10.8.0)
    MESSAGE(FATAL_ERROR "Trilinos version 10.8 or newer required!")
  ENDIF()

  MESSAGE("\nFound Trilinos!  Here are the details: ")
  MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
  MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
  MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
  MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
  MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
  MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
  MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
  MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
  MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
  MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
  MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
  MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
  MESSAGE("End of Trilinos details\n")

  # Get libraries for link line from Trilinos build information
  set(CISM_TRILINOS_LIBS ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Trilinos_EXTRA_LD_FLAGS})

  IF (NOT DEFINED CMAKE_CXX_COMPILER AND NOT DEFINED ENV{CXX}) 
    SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
  ENDIF()
  IF (NOT DEFINED CMAKE_C_COMPILER AND NOT DEFINED ENV{CC}) 
    SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
  ENDIF()
  IF (NOT DEFINED CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC})
    SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
  ENDIF()
ENDIF()

ENABLE_LANGUAGE(Fortran) 

IF (CISM_INCLUDE_IMPLICIT_LINK_LIBRARIES)
  # WJS (6-3-14) Until now, these Fortran_IMPLICIT_LINK_LIBRARIES were always
  # included. However, explicitly appending these implicit link libraries breaks
  # the build on yellowstone, for some reason. It could be because the trilinos
  # build is old. I suspect these implicit link libraries are needed when you're
  # linking with a C++ linker, which is not the case for yellowstone-intel. In
  # any case, I'm providing an option to exclude these from the build, to allow
  # the yellowstone-intel build to work.
  #
  # Also, note that, although these are added to the variable CISM_TRILINOS_LIBS,
  # this variable is used even when building without trilinos. I am merely
  # maintaining the old behavior in this respect.
  LIST(APPEND CISM_TRILINOS_LIBS ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

  #message("")
  #message("  CMake detected the following libraries for linking Fortran with C++ compiler:")
  #message("     ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} ")
ENDIF()

# Only include C++ support if it's really needed, to avoid problems
# caused by broken C++ compilers.
#IK, 8/3/13: added simple_felix option 
IF (CISM_USE_TRILINOS OR CISM_ENABLE_BISICLES)
  ENABLE_LANGUAGE(CXX)
ENDIF()

MESSAGE(">> CISM_NETCDF_DIR set to : ${CISM_NETCDF_DIR}")

IF (NOT DEFINED CISM_NETCDF_DIR)
  MESSAGE(FATAL_ERROR "\nCISM Error: cmake must define CISM_NETCDF_DIR:
  (-D CISM_NETCDF_DIR=<netcdf_install_prefix>)!")
ENDIF()

FIND_PATH(CISM_NETCDFF_FOUND libnetcdff.a ${CISM_NETCDF_DIR}/lib)
#MESSAGE(" CISM_NETCDFF_FOUND = ${CISM_NETCDFF_FOUND}")

IF (${CISM_NETCDFF_FOUND} STREQUAL ${CISM_NETCDF_DIR}/lib )
   SET(CISM_NETCDF_LIBS "netcdff;netcdf" CACHE STRING "Netcdf Library Names(s)")
ELSE()
   SET(CISM_NETCDF_LIBS "netcdf" CACHE STRING "Netcdf Library Names(s)") 
ENDIF()
MESSAGE(">> CISM_NETCDF_LIBS Library(s) set to : ${CISM_NETCDF_LIBS}")

IF (DEFINED CISM_MPI_BASE_DIR)
  IF (NOT DEFINED CISM_MPI_LIB_DIR)
    SET(CISM_MPI_LIB_DIR ${CISM_MPI_BASE_DIR}/lib)
  ENDIF()
  IF (NOT DEFINED CISM_MPI_INC_DIR)
    SET(CISM_MPI_INC_DIR ${CISM_MPI_BASE_DIR}/include)
  ENDIF()
ENDIF()
MESSAGE(">> CISM_MPI_LIB_DIR set to : ${CISM_MPI_LIB_DIR}")
MESSAGE(">> CISM_MPI_INC_DIR set to : ${CISM_MPI_INC_DIR}")


IF (NOT DEFINED CMAKE_Fortran_MODULE_DIRECTORY)
  SET(CMAKE_Fortran_MODULE_DIRECTORY ${CISM_BINARY_DIR}/fortran_mod_files)
ENDIF()
INCLUDE_DIRECTORIES(${CMAKE_Fortran_MODULE_DIRECTORY})


# Note that C++ is excluded here -- we only include C++ support if
# it's really needed (see ENABLE_LANGUAGE(CXX) command above)
PROJECT(CISM Fortran C)


IF (NOT CISM_USE_DEFAULT_IO)
  # Auto-generate the *_io.F90 files in the build directory:
  MESSAGE(">> Calling utils/build/autogenerate-in-build-dir")
  EXECUTE_PROCESS(COMMAND  ${CISM_SOURCE_DIR}/utils/build/autogenerate-in-build-dir
                  ${CISM_SOURCE_DIR} 
                  WORKING_DIRECTORY ${CISM_BINARY_DIR}
                  OUTPUT_FILE ${CISM_BINARY_DIR}/autogenerate.log)
 
  MESSAGE(">> see ${CISM_BINARY_DIR}/autogenerate.log")

  IF (CISM_USE_CISM_FRONT_END)
    # Auto-generate the *_io.F90 files in the build directory:
    MESSAGE(">> Calling utils/build/autogenerate-in-build-dir")
    EXECUTE_PROCESS(COMMAND  ${CISM_SOURCE_DIR}/utils/build/autogen-for-glint-and-glad-in-build-dir
                    ${CISM_SOURCE_DIR} 
                    WORKING_DIRECTORY ${CISM_BINARY_DIR}
                    OUTPUT_FILE ${CISM_BINARY_DIR}/autogen-for-glint-and-glad.log)
 
    MESSAGE(">> see ${CISM_BINARY_DIR}/autogen-for-glint-and-glad.log")
  ENDIF()

ELSE()
  # Simply copy the default io files into the build directory
  MESSAGE(">> Calling utils/build/autocopy-io-to-build-dir")
  EXECUTE_PROCESS(COMMAND  ${CISM_SOURCE_DIR}/utils/build/autocopy-io-to-build-dir
                  ${CISM_SOURCE_DIR}
                  WORKING_DIRECTORY ${CISM_BINARY_DIR}
                  OUTPUT_FILE ${CISM_BINARY_DIR}/autocopy-io.log)
		
  MESSAGE(">> see ${CISM_BINARY_DIR}/autocopy-io.log")
ENDIF()

# Copy a few needed files to the build directory:
MESSAGE(">> Calling utils/build/autocopy-to-build-dir")
EXECUTE_PROCESS(COMMAND  ${CISM_SOURCE_DIR}/utils/build/autocopy-to-build-dir
                ${CISM_SOURCE_DIR}
                WORKING_DIRECTORY ${CISM_BINARY_DIR}
                OUTPUT_FILE ${CISM_BINARY_DIR}/autocopy.log)

MESSAGE(">> see ${CISM_BINARY_DIR}/autocopy.log")
INCLUDE_DIRECTORIES(${CISM_BINARY_DIR}/fortran_autocopy_includes)

# End of setup and error checking
#  NOTE: PROJECT command checks for compilers, so this statement
#        is moved AFTER setting CMAKE_CXX_COMPILER from Trilinos

## Use CMAKE_CXX_FLAGS   CMAKE_Fortran_FLAGS   to override  Trilinos flags
## USe CISM_CXX_FLAGS CISM_Fortran_FLAGS to append to Trilinos flags

IF (NOT CMAKE_CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS} )
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CISM_CXX_FLAGS}")

IF (NOT CMAKE_Fortran_FLAGS)
  SET(CMAKE_Fortran_FLAGS ${Trilinos_Fortran_COMPILER_FLAGS} )
ENDIF()
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CISM_Fortran_FLAGS}")

IF (CISM_USE_TRILINOS) 
   ## CISM requires Trilinos to include the Piro build 
   MESSAGE("-- Looking for Piro within Trilinos:")
   MESSAGE("Trilinos pkg list: " ${Trilinos_PACKAGE_LIST})
   LIST(FIND Trilinos_PACKAGE_LIST  Piro Piro_List_ID)
   IF (Piro_List_ID GREATER -1)
     MESSAGE("-- Looking for Piro: -- found! Continuing.")
   ELSE()
     MESSAGE(FATAL_ERROR "-- Looking for Piro: -- not found! Rebuild Trilinos with Piro package.")
   ENDIF()

   ADD_DEFINITIONS(-DCISM_MPI -DTRILINOS)
    ADD_DEFINITIONS(-DGLIMMER_MPI)
   MESSAGE(">> Trilinos build: adding  -DTRILINOS  to compile lines")

   include_directories(${CISM_NETCDF_DIR}/include)
   link_directories(${CISM_NETCDF_DIR}/lib) 

ELSE()

   IF (${CISM_MPI_MODE})
      ADD_DEFINITIONS(-DCISM_MPI)
      ADD_DEFINITIONS(-DGLIMMER_MPI)
   ENDIF()
   include_directories(${CISM_NETCDF_DIR}/include ${CISM_MPI_INC_DIR})

   link_directories(${CISM_NETCDF_DIR}/lib ${CISM_MPI_LIB_DIR} ${CISM_SCI_LIB_DIR} 
        		 ${CISM_HDF5_LIB_DIR} ) 
ENDIF()

IF (CISM_ENABLE_BISICLES)
  MESSAGE("\nBISICLES build configured.\n")  

  SET(CISM_DYCORE_DIR ${CISM_SOURCE_DIR}/libdycore)
  SET(CISM_TO_DYCORE ${CISM_DYCORE_DIR}/glimmer_to_dycore.F90) 
  SET(CISM_BISICLES_DIR ${CISM_DYCORE_DIR}/BISICLES)
  
  INCLUDE_DIRECTORIES(${CISM_DYCORE_DIR})
  LINK_DIRECTORIES(${CISM_DYCORE_DIR} ${BISICLES_INTERFACE_DIR} )
  ADD_DEFINITIONS(-DCISM_HAS_BISICLES)
ENDIF()


IF ( CISM_USE_GPTL_INSTRUMENTATION AND 
    (NOT CISM_SERIAL_MODE))
 OPTION(CISM_GPTL_DIR "Optional path to installed gptl library")
 IF (CISM_GPTL_DIR)
   message(">> GPTL Library reqested, installed in ${CISM_GPTL_DIR}")
   ADD_DEFINITIONS(-DCESMTIMERS)
   include_directories(${CISM_GPTL_DIR})
   link_directories(${CISM_GPTL_DIR})
   SET(CISM_GPTL_LIB gptl)
 ELSE()
   message(">>GPTL Library not requested: can set CISM_GPTL_DIR to enable")
 ENDIF()
ENDIF()

#####Build All Fortran Sources #####

#IF (CISM_USE_CISM_FRONT_END)
   file(GLOB GLINT_SOURCES libglint/*.F90)
   file(GLOB GLAD_SOURCES libglad/*.F90)
   file(GLOB GLISSADE_SOURCES libglissade/*.F90)
#ENDIF()


file(GLOB FORTRANSOURCES
     libglimmer-solve/SLAP/*.f
     libglimmer-solve/*.F90
     libglimmer/*.F90 libglimmer/writestats.c
     libglide/*.F90
     cism_driver/eismint_forcing.F90
     cism_driver/testsfg.F90
     ${GLINT_SOURCES}
     ${GLAD_SOURCES}
     ${GLISSADE_SOURCES}
     ${CISM_TO_DYCORE})

LIST(REMOVE_ITEM FORTRANSOURCES
   ${CISM_SOURCE_DIR}/libglimmer-solve/SLAP/dlapqc.f )


IF (CISM_USE_TRILINOS OR CISM_MPI_MODE)
   LIST(REMOVE_ITEM FORTRANSOURCES
      ${CISM_SOURCE_DIR}/libglimmer/parallel_slap.F90)   

ELSE(${CISM_SERIAL_MODE})
   LIST(REMOVE_ITEM FORTRANSOURCES
      ${CISM_SOURCE_DIR}/libglimmer/parallel_mpi.F90)

   IF (CISM_USE_MPI_WITH_SLAP)
     ADD_DEFINITIONS(-D_USE_MPI_WITH_SLAP)
   ELSE()
     LIST(REMOVE_ITEM FORTRANSOURCES
       ${CISM_SOURCE_DIR}/libglimmer/mpi_mod.F90)
   ENDIF()
ENDIF()

# Get autogenerated source files, and add them to the FORTRANSOURCES list:
# changed to always do this:
#IF (CISM_COUPLED)
  FILE(GLOB FORTRAN_AUTOGEN_SOURCES ${CISM_BINARY_DIR}/fortran_autogen_srcs/*.F90) 
  message("Autogenerated CISM sources: ${FORTRAN_AUTOGEN_SOURCES}")

  SET(FORTRANSOURCES ${FORTRANSOURCES} ${FORTRAN_AUTOGEN_SOURCES})

# Remove old versions of autogenerated F90 files that may be sitting
# around in the source tree from an old cmake-based build or an
# autotools-based build. (Now the cmake-based build doesn't place
# these in the source tree.)
# Unlike SOURCEMODFILES, we just hard-code the files to remove for
# simplicity.
# This can be removed once we switch to consistently using this new
# cmake-based build.
  LIST(REMOVE_ITEM FORTRANSOURCES
     ${CISM_SOURCE_DIR}/libglide/glide_io.F90
     ${CISM_SOURCE_DIR}/libglide/glide_lithot_io.F90
     ${CISM_SOURCE_DIR}/libglint/glint_io.F90
     ${CISM_SOURCE_DIR}/libglint/glint_mbal_io.F90
     ${CISM_SOURCE_DIR}/libglimmer/glimmer_vers.F90 )
#ENDIF()


### (For CESM) Remove source files with names already in CISM_SOURCEMOD_DIR
OPTION(CISM_SOURCEMOD_DIR 
      "Path to SourceMod directory of F90 files to replace Glimmer files")

#MESSAGE("Fortran Source Files: ${FORTRANSOURCES}")

# Note that the following glob does NOT contain .cpp files, because
# those are built in a separate library - so for now, you can't put
# .cpp files in your sourceMod directory.
FILE(GLOB SOURCEMODFILES 
     ${CISM_SOURCEMOD_DIR}/*.F90
     ${CISM_SOURCEMOD_DIR}/*.F
     ${CISM_SOURCEMOD_DIR}/*.f90
     ${CISM_SOURCEMOD_DIR}/*.f
     ${CISM_SOURCEMOD_DIR}/*.c)

# MESSAGE("glimmer_sourcemod_dir: " ${CISM_SOURCEMOD_DIR})
# MESSAGE("Fortran Mod Files: ${SOURCEMODFILES}")

FOREACH( MODFILE ${SOURCEMODFILES})
  STRING(FIND ${MODFILE} / index REVERSE)
  MATH(EXPR index ${index}+1)
  STRING(SUBSTRING ${MODFILE} ${index} -1 filename)
  FOREACH( SOURCEFILE ${FORTRANSOURCES})
    STRING(REGEX MATCH ${filename} match_found ${SOURCEFILE})

    IF(match_found)
      MESSAGE("--SourceMod: removing ${SOURCEFILE} in favor of ${MODFILE}")
      LIST(REMOVE_ITEM FORTRANSOURCES ${SOURCEFILE})
   ENDIF()
  ENDFOREACH()
ENDFOREACH()
###

#IF (CISM_COUPLED)
  # enable removal of the autogenerated source files, when 'make clean' is done:
  # commented out for now, since it may create more problems than it solves when building
  # SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${FORTRAN_AUTOGEN_SOURCES}") 
#ENDIF()


SET(CMAKE_INCLUDE_PATH ${CISM_BINARY_DIR}/include)


#message("Fortransource ${FORTRANSOURCES}")

# include-dirs: Binary: all the .mod files; SOURCE: config.inc; 
#               libglide: glide_mask.inc

include_directories (${CISM_BINARY_DIR} ${CISM_SOURCE_DIR}
                     ${CISM_SOURCE_DIR}/libglide)


SET(LIBRARY_OUTPUT_PATH ${CISM_BINARY_DIR}/lib)

### Determine which fortran source files are fixed-form,
### which may require special compilation flags
FOREACH (SOURCEFILE ${FORTRANSOURCES} ${SOURCEMODFILES})
  STRING(REGEX MATCH "\\.[fF]$" match_found ${SOURCEFILE})
  IF(match_found)
    LIST(APPEND FIXEDSOURCES ${SOURCEFILE})
  ENDIF()
ENDFOREACH()
SET_SOURCE_FILES_PROPERTIES(${FIXEDSOURCES} PROPERTIES Fortran_FORMAT FIXED)

add_library(glimmercismfortran ${FORTRANSOURCES} ${SOURCEMODFILES})


#####Build C++ Sources #####

IF (CISM_USE_TRILINOS)
  add_subdirectory(libglimmer-trilinos)
ENDIF()

#####Build simple_glide executable (unless turned off) #####
OPTION(CISM_NO_EXECUTABLE "Set to  ON  to just build libraries (default:OFF)"  OFF)


#IF (CISM_BUILD_SIMPLE_GLIDE AND (NOT CISM_NO_EXECUTABLE))
#   add_subdirectory(example-drivers/simple_glide/src)
#ENDIF()

IF (CISM_USE_CISM_FRONT_END)
   add_subdirectory(cism_driver)
ENDIF()
 

IF (CISM_ENABLE_BISICLES)
  add_subdirectory(libdycore)

  include_directories (${CISM_DYCORE_DIR}
                       ${CISM_BISICLES_DIR})

  message("glimmer src dir: ${CISM_SOURCE_DIR}")
  message("glimmer dycore dir: ${CISM_DYCORE_DIR}")
  get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
  message("inc_dirs = ${inc_dirs}")
ENDIF()


