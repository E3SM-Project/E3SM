# CMAKE File for GLIMMER building against an installed Trilinos

cmake_minimum_required(VERSION 2.8.6)

OPTION (NO_TRILINOS "Toggle to Not use Trilinos: defaults to OFF (Trilinos on by default)" OFF)
OPTION (CISM_MPI_MODE "Toggle to Configure with MPI: defaults to ON" ON)
OPTION (CISM_SERIAL_MODE "Toggle to Configure in Serial mode: defaults to OFF " OFF)
OPTION (CISM_USE_MPI_WITH_SLAP "Toggle to use mpi when using SLAP solver, only relevant if CISM_SERIAL_MODE=ON: defaults to OFF" OFF)
OPTION (CISM_BUILD_SIMPLE_GLIDE "Toggle to build simple_glide, ON by default" ON)
OPTION (CISM_BUILD_SIMPLE_BISICLES "Toggle to build simple_bisicles, off by default" OFF)
OPTION (CISM_BUILD_EXTRA_EXECUTABLES "Toggle to other executables, off by default" OFF)
OPTION (CISM_USE_GPTL_INSTRUMENTATION "Toggle to use GPTL instrumentation, on by default " ON)
OPTION (CISM_BUILD_GLINT "Toggle to build glint, off by default" OFF)
OPTION (CISM_COUPLED "Toggle to build CISM for use with CESM, off by default" OFF)
OPTION (CISM_USE_CISM_FRONT_END "Toggle to use cism_driver or cism_cesm_interface with cism_front_end, off by default" OFF)
OPTION (CISM_USE_DEFAULT_IO "Toggle to use default i/o files rather than running python script, off by default" OFF)
# WJS (1-3-13): We could avoid CISM_GNU by using something like CMAKE_Fortran_COMPILER_ID or CMAKE_COMPILER_IS_GNUCC, 
# but it's not clear to me if those work consistently
OPTION (CISM_GNU "Toggle to set compilation flags needed for the gnu compiler, off by default" OFF)

# override CISM_SERIAL_MODE setting, if it conflicts with CISM_MPI_MODE setting:
IF (CISM_MPI_MODE)
  SET(CISM_SERIAL_MODE "OFF")
  SET(CISM_BUILD_EXTRA_EXECUTABLES "OFF")
ENDIF()

IF (NOT CISM_MPI_MODE AND NOT CISM_SERIAL_MODE)
  SET(CISM_MPI_MODE "ON")
ENDIF()

#MESSAGE("CISM MPI Mode Setting:    ${CISM_MPI_MODE}")
#MESSAGE("CISM Serial Mode Setting: ${CISM_SERIAL_MODE}")


# set (or override) other options, if CISM_COUPLED is ON:
IF (CISM_COUPLED)
  SET(GLIMMER_NO_EXECUTABLE "ON")
  SET(CISM_BUILD_SIMPLE_GLIDE "OFF")
  SET(CISM_BUILD_SIMPLE_BISICLES "OFF")
  SET(CISM_BUILD_EXTRA_EXECUTABLES "OFF")
  SET(CISM_BUILD_GLINT "ON")  
  SET(CISM_USE_DEFAULT_IO "ON")
  ADD_DEFINITIONS(-DCCSMCOUPLED)
ENDIF()


IF (NOT DEFINED CISM_BINARY_DIR)
   SET(CISM_BINARY_DIR ${CMAKE_BINARY_DIR})
ENDIF()

MESSAGE("Building in: ${CISM_BINARY_DIR}")

IF (NOT ${NO_TRILINOS})
  OPTION(GLIMMER_TRILINOS_DIR "Required path to installed Trilinos")
  OPTION(GLIMMER_NETCDF_DIR "Required path to installed Netcdf")

  # Error check up front
  IF (NOT DEFINED GLIMMER_TRILINOS_DIR)
    MESSAGE(FATAL_ERROR "\nGLIMMER Error: cmake must define GLIMMER_TRILINOS_DIR:
    (-D GLIMMER_TRILINOS_DIR=<trilinos_install_prefix>)!")
  ENDIF()

  # Get Trilinos as one entity
  SET(CMAKE_PREFIX_PATH ${GLIMMER_TRILINOS_DIR} ${CMAKE_PREFIX_PATH})
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
# Only include C++ support if it's really needed, to avoid problems
# caused by broken C++ compilers.
IF ((NOT NO_TRILINOS) OR CISM_BUILD_SIMPLE_BISICLES)
  ENABLE_LANGUAGE(CXX)
ENDIF()

MESSAGE(">> GLIMMER_NETCDF_DIR set to : ${GLIMMER_NETCDF_DIR}")

IF (NOT DEFINED GLIMMER_NETCDF_DIR)
  MESSAGE(FATAL_ERROR "\nGLIMMER Error: cmake must define GLIMMER_NETCDF_DIR:
  (-D GLIMMER_NETCDF_DIR=<netcdf_install_prefix>)!")
ENDIF()


FIND_PATH(CISM_NETCDFF_FOUND libnetcdff.a ${GLIMMER_NETCDF_DIR}/lib)

#MESSAGE(" CISM_NETCDFF_FOUND = ${CISM_NETCDFF_FOUND}")
IF (${CISM_NETCDFF_FOUND} STREQUAL ${GLIMMER_NETCDF_DIR}/lib )
   SET(GLIMMER_NETCDF_LIBS "netcdff;netcdf" CACHE STRING "Netcdf Library Names(s)")
ELSE()
   SET(GLIMMER_NETCDF_LIBS "netcdf" CACHE STRING "Netcdf Library Names(s)") 
ENDIF()
MESSAGE(">> GLIMMER_NETCDF_LIBS Library(s) set to : ${GLIMMER_NETCDF_LIBS}")

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
PROJECT(GLIMMER Fortran C)


IF (NOT CISM_USE_DEFAULT_IO)
  # Auto-generate the *_io.F90 files in the build directory:
  MESSAGE(">> Calling cmake-scripts/autogenerate-in-build-dir")
  EXECUTE_PROCESS(COMMAND  ${GLIMMER_SOURCE_DIR}/cmake-scripts/autogenerate-in-build-dir
                  ${GLIMMER_SOURCE_DIR} 
                  WORKING_DIRECTORY ${CISM_BINARY_DIR}
                  OUTPUT_FILE ${CISM_BINARY_DIR}/autogenerate.log)
  MESSAGE(">> see ${CISM_BINARY_DIR}/autogenerate.log")
ELSE()
  # Simply copy the default io files into the build directory
  MESSAGE(">> Calling cmake-scripts/autocopy-io-to-build-dir")
  EXECUTE_PROCESS(COMMAND  ${GLIMMER_SOURCE_DIR}/cmake-scripts/autocopy-io-to-build-dir
                  ${GLIMMER_SOURCE_DIR}
                  WORKING_DIRECTORY ${CISM_BINARY_DIR}
                  OUTPUT_FILE ${CISM_BINARY_DIR}/autocopy-io.log)
		
  MESSAGE(">> see ${CISM_BINARY_DIR}/autocopy-io.log")
ENDIF()

# Copy a few needed files to the build directory:
MESSAGE(">> Calling cmake-scripts/autocopy-to-build-dir")
EXECUTE_PROCESS(COMMAND  ${GLIMMER_SOURCE_DIR}/cmake-scripts/autocopy-to-build-dir
                ${GLIMMER_SOURCE_DIR}
                WORKING_DIRECTORY ${CISM_BINARY_DIR}
                OUTPUT_FILE ${CISM_BINARY_DIR}/autocopy.log)

MESSAGE(">> see ${CISM_BINARY_DIR}/autocopy.log")
INCLUDE_DIRECTORIES(${CISM_BINARY_DIR}/fortran_autocopy_includes)

# End of setup and error checking
#  NOTE: PROJECT command checks for compilers, so this statement
#        is moved AFTER setting CMAKE_CXX_COMPILER from Trilinos

## Use CMAKE_CXX_FLAGS   CMAKE_Fortran_FLAGS   to override  Trilinos flags
## USe GLIMMER_CXX_FLAGS GLIMMER_Fortran_FLAGS to append to Trilinos flags

IF (NOT CMAKE_CXX_FLAGS)
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS} )
ENDIF()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GLIMMER_CXX_FLAGS}")

IF (NOT CMAKE_Fortran_FLAGS)
  SET(CMAKE_Fortran_FLAGS ${Trilinos_Fortran_COMPILER_FLAGS} )
ENDIF()
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GLIMMER_Fortran_FLAGS}")

IF (NOT ${NO_TRILINOS}) 
   ## GLIMMER requires Trilinos to include the Piro build 
   MESSAGE("-- Looking for Piro within Trilinos:")
   MESSAGE("Trilinos pkg list: " ${Trilinos_PACKAGE_LIST})
   LIST(FIND Trilinos_PACKAGE_LIST  Piro Piro_List_ID)
   IF (Piro_List_ID GREATER -1)
     MESSAGE("-- Looking for Piro: -- found! Continuing.")
   ELSE()
     MESSAGE(FATAL_ERROR "-- Looking for Piro: -- not found! Rebuild Trilinos with Piro package.")
   ENDIF()

   ADD_DEFINITIONS(-DGLIMMER_MPI -DTRILINOS)
   MESSAGE(">> Trilinos build: adding  -DTRILINOS  to compile lines")

   include_directories(${GLIMMER_NETCDF_DIR}/include)
   link_directories(${GLIMMER_NETCDF_DIR}/lib) 

ELSE()

   IF (${CISM_MPI_MODE})
      ADD_DEFINITIONS(-DGLIMMER_MPI)
   ENDIF()
   include_directories(${GLIMMER_NETCDF_DIR}/include ${CISM_MPI_INC_DIR})

#   link_directories(${GLIMMER_NETCDF_DIR}/lib )
   link_directories(${GLIMMER_NETCDF_DIR}/lib ${CISM_MPI_LIB_DIR} ${CISM_SCI_LIB_DIR} 
        		 ${CISM_HDF5_LIB_DIR} ) 

ENDIF()

IF (CISM_BUILD_SIMPLE_BISICLES)
  MESSAGE("\nBISICLES build configured.\n")  

  SET(GLIMMER_DYCORE_DIR ${GLIMMER_SOURCE_DIR}/libdycore)
  SET(GLIMMER_TO_DYCORE ${GLIMMER_DYCORE_DIR}/glimmer_to_dycore.F90) 
  SET(GLIMMER_BISICLES_DIR ${GLIMMER_DYCORE_DIR}/BISICLES)
  
  INCLUDE_DIRECTORIES(${GLIMMER_DYCORE_DIR})
  LINK_DIRECTORIES(${GLIMMER_DYCORE_DIR} ${BISICLES_INTERFACE_DIR} )
ENDIF()

IF ( CISM_USE_GPTL_INSTRUMENTATION AND 
    (NOT CISM_BUILD_SIMPLE_BISICLES) AND
    (NOT CISM_SERIAL_MODE))
 OPTION(GLIMMER_GPTL_DIR "Optional path to installed gptl library")
 IF (GLIMMER_GPTL_DIR)
   message(">> GPTL Library reqested, installed in ${GLIMMER_GPTL_DIR}")
   ADD_DEFINITIONS(-DCESMTIMERS)
   include_directories(${GLIMMER_GPTL_DIR})
   link_directories(${GLIMMER_GPTL_DIR})
   SET(GLIMMER_GPTL_LIB gptl)
 ELSE()
   message(">>GPTL Library not requested: can set GLIMMER_GPTL_DIR to enable")
 ENDIF()
ENDIF()

#####Build All Fortran Sources #####

IF (CISM_BUILD_GLINT OR CISM_BUILD_EXTRA_EXECUTABLES)
   file(GLOB GLINT_SOURCES libglint/*.F90) 
ENDIF()



file(GLOB FORTRANSOURCES
     libglimmer-solve/SLAP/*.f
     libglimmer-solve/*.F90
     libglimmer/*.F90 libglimmer/writestats.c
     libglide/*.F90
     ${GLINT_SOURCES}
     ${GLIMMER_TO_DYCORE})

LIST(REMOVE_ITEM FORTRANSOURCES
   ${GLIMMER_SOURCE_DIR}/libglimmer/unittest_glimmer_deriv.F90 
   ${GLIMMER_SOURCE_DIR}/libglide/relaxed.F90
   ${GLIMMER_SOURCE_DIR}/libglide/test_setup.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_commandline.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_config.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_filenames.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_integrate.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_ts.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/test_writestats.F90
   ${GLIMMER_SOURCE_DIR}/libglimmer/nc2config.F90 )


IF ((NOT ${NO_TRILINOS}) OR CISM_MPI_MODE)
   LIST(REMOVE_ITEM FORTRANSOURCES
      ${GLIMMER_SOURCE_DIR}/libglimmer/parallel_slap.F90   
      ${GLIMMER_SOURCE_DIR}/libglimmer/glimmer_horiz_bcs_serial.F90)

ELSE(${CISM_SERIAL_MODE})
   LIST(REMOVE_ITEM FORTRANSOURCES
      ${GLIMMER_SOURCE_DIR}/libglimmer/parallel_mpi.F90
      ${GLIMMER_SOURCE_DIR}/libglimmer/glimmer_horiz_bcs_parallel.F90)

   IF (CISM_USE_MPI_WITH_SLAP)
     ADD_DEFINITIONS(-D_USE_MPI_WITH_SLAP)
   ELSE()
     LIST(REMOVE_ITEM FORTRANSOURCES
       ${GLIMMER_SOURCE_DIR}/libglimmer/mpi_mod.F90)
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
     ${GLIMMER_SOURCE_DIR}/libglide/glide_io.F90
     ${GLIMMER_SOURCE_DIR}/libglide/glide_lithot_io.F90
     ${GLIMMER_SOURCE_DIR}/libglint/glint_io.F90
     ${GLIMMER_SOURCE_DIR}/libglint/glint_mbal_io.F90
     ${GLIMMER_SOURCE_DIR}/libglimmer/glimmer_vers.F90 
     ${GLIMMER_SOURCE_DIR}/example_drivers/eis/src/eis_io.F90 )
#ENDIF()


### (For CESM) Remove source files with names already in GLIMMER_SOURCEMOD_DIR
OPTION(GLIMMER_SOURCEMOD_DIR 
      "Path to SourceMod directory of F90 files to replace Glimmer files")

#MESSAGE("Fortran Source Files: ${FORTRANSOURCES}")

# Note that the following glob does NOT contain .cpp files, because
# those are built in a separate library - so for now, you can't put
# .cpp files in your sourceMod directory.
FILE(GLOB SOURCEMODFILES 
     ${GLIMMER_SOURCEMOD_DIR}/*.F90
     ${GLIMMER_SOURCEMOD_DIR}/*.F
     ${GLIMMER_SOURCEMOD_DIR}/*.f90
     ${GLIMMER_SOURCEMOD_DIR}/*.f
     ${GLIMMER_SOURCEMOD_DIR}/*.c)

# MESSAGE("glimmer_sourcemod_dir: " ${GLIMMER_SOURCEMOD_DIR})
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

# WJS (1-3-13): Ideally, rather than checking CISM_GNU, we would instead check for whether -fno-range-check works,
# and if so, include that flag (see 
# http://stackoverflow.com/questions/3134660/how-to-apply-different-compiler-options-for-different-compilers-in-cmake).
# But it doesn't look like there is that capability for fortran compilers yet.
IF (CISM_GNU)
  # Allow explicit NaN values in gfortran compiler
  # Note that this won't work if the user has put nan_mod.F90 in their
  # SourceMods directory, since it assumes a particular path
  SET_PROPERTY(SOURCE ${GLIMMER_SOURCE_DIR}/libglimmer/nan_mod.F90 
               APPEND PROPERTY COMPILE_FLAGS -fno-range-check)
ENDIF()

#IF (CISM_COUPLED)
  # enable removal of the autogenerated source files, when 'make clean' is done:
  # commented out for now, since it may create more problems than it solves when building
  # SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${FORTRAN_AUTOGEN_SOURCES}") 
#ENDIF()


SET(CMAKE_INCLUDE_PATH ${CISM_BINARY_DIR}/include)


#message("Fortransource ${FORTRANSOURCES}")

# include-dirs: Binary: all the .mod files; SOURCE: config.inc; 
#               libglide: glide_mask.inc

include_directories (${CISM_BINARY_DIR} ${GLIMMER_SOURCE_DIR}
                     ${GLIMMER_SOURCE_DIR}/libglide)


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

IF (NOT ${NO_TRILINOS})
  add_subdirectory(libglimmer-trilinos)
ENDIF()

#####Build simple_glide executable (unless turned off) #####
OPTION(GLIMMER_NO_EXECUTABLE "Set to  ON  to just build libraries (default:OFF)"  OFF)


IF (CISM_BUILD_SIMPLE_GLIDE AND (NOT GLIMMER_NO_EXECUTABLE))
  add_subdirectory(example-drivers/simple_glide/src)
ENDIF()


IF (CISM_BUILD_SIMPLE_BISICLES)
  add_subdirectory(example-drivers/simple_bisicles/src)
  add_subdirectory(libdycore)

  IF (CISM_USE_CISM_FRONT_END)
    add_subdirectory(cism_driver)
  ENDIF()

  SET(GLIMMER_EXTERNAL_DYCORE_DIR ${GLIMMER_SOURCE_DIR}/libdycore)
  SET(GLIMMER_BISICLES_DIR ${GLIMMER_EXTERNAL_DYCORE_DIR}/BISICLES)
  include_directories (${GLIMMER_EXTERNAL_DYCORE_DIR}
                       ${GLIMMER_BISICLES_DIR})

  message("glimmer src dir: ${GLIMMER_SOURCE_DIR}")
  message("glimmer extern dycore dir: ${GLIMMER_EXTERNAL_DYCORE_DIR}")
  get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
  message("inc_dirs = ${inc_dirs}")

ENDIF()


IF(CISM_BUILD_EXTRA_EXECUTABLES AND CISM_SERIAL_MODE)
  add_subdirectory(example-drivers/eis/src)
  add_subdirectory(example-drivers/glint/src)
  message("Building extra executables: ")
  MESSAGE(">> Calling cmake-scripts/autogen-extras-in-build-dir")
  EXECUTE_PROCESS(COMMAND  ${GLIMMER_SOURCE_DIR}/cmake-scripts/autogen-extras-in-build-dir
                  ${GLIMMER_SOURCE_DIR} 
                  WORKING_DIRECTORY ${CISM_BINARY_DIR}
                  OUTPUT_FILE ${CISM_BINARY_DIR}/autogen-extras.log)
  MESSAGE(">> see ${CISM_BINARY_DIR}/autogen-extras.log")
ENDIF()

#COMMAND (environment)


