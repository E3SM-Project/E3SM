function (prc var)
  message("${var} ${${var}}")
endfunction ()

# Macro to create config file. The macro creates a temporary config
# file first. If the config file in the build directory does not
# exist or it is different from the temporary one, then the config
# file is updated. This avoid rebuilding the (almost) entire homme
# every time cmake is run.

MACRO (HommeConfigFile CONFIG_FILE_IN CONFIG_FILE_C CONFIG_FILE_F90)

  CONFIGURE_FILE (${CONFIG_FILE_IN} ${CONFIG_FILE_C}.tmp)

  # Assume by default that config file is out of date
  SET (OUT_OF_DATE TRUE)

  # If config file in binary dir exists, we check whether the new one would be different
  IF (EXISTS ${CONFIG_FILE_C})

    # We rely on FILE macro rather than running diff, since it is
    # more portable (guaranteed to work regardless of underlying system)
    FILE (READ ${CONFIG_FILE_C} CONFIG_FILE_C_STR)
    FILE (READ ${CONFIG_FILE_C}.tmp CONFIG_FILE_C_TMP_STR)

    IF (CONFIG_FILE_C_STR STREQUAL CONFIG_FILE_C_TMP_STR)
      # config file was present and appears unchanged
      SET (OUT_OF_DATE FALSE)
    ENDIF()

    FILE (REMOVE ${CONFIG_FILE_C}.tmp)
  ENDIF ()

  # If out of date (either missing or different), adjust
  IF (OUT_OF_DATE)

    # Run the configure macro
    CONFIGURE_FILE (${CONFIG_FILE_IN} ${CONFIG_FILE_C})

    # Run sed to change '/*...*/' comments into '!/*...*/'
    EXECUTE_PROCESS(COMMAND sed "s;^/;!/;g"
                    WORKING_DIRECTORY ${HOMME_BINARY_DIR}
                    INPUT_FILE ${CONFIG_FILE_C}
                    OUTPUT_FILE ${CONFIG_FILE_F90})
  ENDIF()

ENDMACRO (HommeConfigFile)

# Macro to create the individual tests
macro(createTestExec execName execType macroNP macroNC
                     macroPLEV macroUSE_PIO macroWITH_ENERGY macroQSIZE_D)

# before calling this macro, be sure that these are set locally:
# EXEC_INCLUDE_DIRS
# EXEC_SOURCES

  # Set the include directories
  SET(EXEC_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${execName}_modules")
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}
                      ${EXEC_INCLUDE_DIRS}
                      ${EXEC_MODULE_DIR}
                      )

  MESSAGE(STATUS "Building ${execName} derived from ${execType} with:")
  MESSAGE(STATUS "  NP = ${macroNP}")
  MESSAGE(STATUS "  PLEV = ${macroPLEV}")
  MESSAGE(STATUS "  QSIZE_D = ${macroQSIZE_D}")
  MESSAGE(STATUS "  PIO = ${macroUSE_PIO}")
  MESSAGE(STATUS "  ENERGY = ${macroWITH_ENERGY}")

  # Set the variable to the macro variables
  SET(NUM_POINTS ${macroNP})
  SET(NUM_CELLS ${macroNC})
  SET(NUM_PLEV ${macroPLEV})

  IF (${macroUSE_PIO})
    SET(PIO TRUE)
    SET(PIO_INTERP)
  ELSE ()
    SET(PIO)
    SET(PIO_INTERP TRUE)
  ENDIF ()

  IF(BUILD_HOMME_WITHOUT_PIOLIBRARY AND (NOT PIO_INTERP))
    MESSAGE(ERROR "For building without PIO library set PREQX_USE_PIO to false")
  ENDIF ()

  IF (${macroWITH_ENERGY})
    SET(ENERGY_DIAGNOSTICS TRUE)
  ELSE()
    SET(ENERGY_DIAGNOSTICS)
  ENDIF ()

  IF (${macroQSIZE_D})
    SET(QSIZE_D ${macroQSIZE_D})
  ELSE()
    SET(QSIZE_D)
  ENDIF ()


  # This is needed to compile the test executables with the correct options
  SET(THIS_CONFIG_IN ${HOMME_SOURCE_DIR}/src/${execType}/config.h.cmake.in)
  SET(THIS_CONFIG_HC ${CMAKE_CURRENT_BINARY_DIR}/config.h.c)
  SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)

  # First configure the file (which formats the file as C)
  HommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

  ADD_DEFINITIONS(-DHAVE_CONFIG_H)

  ADD_EXECUTABLE(${execName} ${EXEC_SOURCES})
  # For SYCL builds it is suggested to use CXX linker with `-fortlib`
  # for mixed-language setups
  IF(Kokkos_ENABLE_SYCL)
    SET_TARGET_PROPERTIES(${execName} PROPERTIES LINKER_LANGUAGE CXX)
  ELSE()
    SET_TARGET_PROPERTIES(${execName} PROPERTIES LINKER_LANGUAGE Fortran)
  ENDIF()
  IF(BUILD_HOMME_WITHOUT_PIOLIBRARY)
    TARGET_COMPILE_DEFINITIONS(${execName} PUBLIC HOMME_WITHOUT_PIOLIBRARY)
  ENDIF()
  IF(BUILD_HOMMEXX_BENCHMARK_NOFORCING)
    TARGET_COMPILE_DEFINITIONS(${execName} PUBLIC HOMMEXX_BENCHMARK_NOFORCING)
  ENDIF()

  target_link_libraries(${execName} csm_share)

  IF (CXXLIB_SUPPORTED_CACHE)
    MESSAGE(STATUS "   Linking Fortran with -cxxlib")
    TARGET_LINK_LIBRARIES(${execName} -cxxlib)
  ENDIF ()

  STRING(TOUPPER "${PERFORMANCE_PROFILE}" PERF_PROF_UPPER)
  IF ("${PERF_PROF_UPPER}" STREQUAL "VTUNE")
    TARGET_LINK_LIBRARIES(${execName} ittnotify)
  ENDIF ()

  # Add this executable to a list
  SET(EXEC_LIST ${EXEC_LIST} ${execName} CACHE INTERNAL "List of configured executables")

  # Record which model this executable runs (e.g., theta-l or theta-l_kokkos)
  SET_TARGET_PROPERTIES(${execName} PROPERTIES HOMME_MODEL ${execType})

  # If this is a Kokkos executable, e.g. theta-l_kokkos, then link to the C++
  # Compose library; if not, then link to the F90 one.
  #   If Compose is not enabled, then COMPOSE_LIBRARY_F90 and
  # COMPOSE_LIBRARY_CPP are empty, so then COMPOSE_LIBRARY_TYPE will be, too.
  string(FIND ${execType} "kokkos" KOKKOS_SUFFIX_LOC)
  if (KOKKOS_SUFFIX_LOC EQUAL -1)
    set (COMPOSE_LIBRARY_TYPE ${COMPOSE_LIBRARY_F90})
  else ()
    set (COMPOSE_LIBRARY_TYPE ${COMPOSE_LIBRARY_CPP})
  endif ()

  TARGET_LINK_LIBRARIES(${execName} timing ${COMPOSE_LIBRARY_TYPE} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
  IF(NOT BUILD_HOMME_WITHOUT_PIOLIBRARY)
    IF(HOMME_USE_SCORPIO)
      TARGET_LINK_LIBRARIES(${execName} piof pioc)
    ELSE ()
      TARGET_LINK_LIBRARIES(${execName} pio)
    ENDIF ()
  ENDIF ()

  IF (HOMME_USE_KOKKOS)
    TARGET_LINK_LIBRARIES(${execName} Kokkos::kokkos)
  ENDIF ()

  # Move the module files out of the way so the parallel build
  # doesn't have a race condition
  SET_TARGET_PROPERTIES(${execName}
    PROPERTIES Fortran_MODULE_DIRECTORY ${EXEC_MODULE_DIR})

  IF (HOMME_USE_MKL)
    IF (MKL_TYPE STREQUAL "oneMKL")
      TARGET_LINK_LIBRARIES(${execName} -qmkl)
    ELSEIF (MKL_TYPE STREQUAL "Intel MKL")
      TARGET_LINK_LIBRARIES(${execName} -mkl)
    ENDIF ()
  ELSE()
    IF (NOT HOMME_FIND_BLASLAPACK)
      TARGET_LINK_LIBRARIES(${execName} lapack blas)
      ADD_DEPENDENCIES(${execName} blas lapack)
    ENDIF()
  ENDIF()

  IF (HAVE_EXTRAE)
    TARGET_LINK_LIBRARIES(${execName} ${Extrae_LIBRARY})
  ENDIF ()

  IF (HOMME_USE_TRILINOS)
    TARGET_LINK_LIBRARIES(${execName} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  ENDIF()


  IF (HOMME_USE_ARKODE AND "${execType}" STREQUAL "theta-l")
    TARGET_LINK_LIBRARIES(${execName} sundials_farkode_mod)
    TARGET_LINK_LIBRARIES(${execName} sundials_arkode)
  ENDIF ()

  INSTALL(TARGETS ${execName} RUNTIME DESTINATION tests)

endmacro(createTestExec)

# Create a library instead of an executable, so Homme can be used by
# another cmake project as a dycore library
macro(createExecLib libName execType libSrcs inclDirs macroNP
                    macroPLEV macroWITH_ENERGY macroQSIZE_D)

  # Set the include directories
  SET(modulesDir "${CMAKE_CURRENT_BINARY_DIR}/${libName}_modules")

  MESSAGE(STATUS "Building ${libName} library derived from ${execType} with:")
  MESSAGE(STATUS "  NP = ${macroNP}")
  MESSAGE(STATUS "  PLEV = ${macroPLEV}")
  MESSAGE(STATUS "  QSIZE_D = ${macroQSIZE_D}")
  MESSAGE(STATUS "  ENERGY = ${macroWITH_ENERGY}")

  # Set the variable to the macro variables
  SET(NUM_POINTS ${macroNP})
  SET(NUM_PLEV ${macroPLEV})

  IF (${macroWITH_ENERGY})
    SET(ENERGY_DIAGNOSTICS TRUE)
  ELSE()
    SET(ENERGY_DIAGNOSTICS)
  ENDIF ()

  IF (${macroQSIZE_D})
    SET(QSIZE_D ${macroQSIZE_D})
  ELSE()
    SET(QSIZE_D)
  ENDIF ()

  # This is needed to compile the test executables with the correct options
  SET(THIS_CONFIG_IN ${HOMME_SOURCE_DIR}/src/${execType}/config.h.cmake.in)
  SET(THIS_CONFIG_HC ${CMAKE_CURRENT_BINARY_DIR}/config.h.c)
  SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)

  # First configure the file (which formats the file as C)
  HommeConfigFile (${THIS_CONFIG_IN} ${THIS_CONFIG_HC} ${THIS_CONFIG_H} )

  ADD_DEFINITIONS(-DHAVE_CONFIG_H)

  ADD_LIBRARY(${libName} ${libSrcs})
  TARGET_INCLUDE_DIRECTORIES (${libName} PUBLIC ${inclDirs} ${modulesDir} ${CMAKE_CURRENT_BINARY_DIR})
  SET_TARGET_PROPERTIES(${libName} PROPERTIES Fortran_MODULE_DIRECTORY ${modulesDir})
  SET_TARGET_PROPERTIES(${libName} PROPERTIES LINKER_LANGUAGE Fortran)
  IF(BUILD_HOMME_WITHOUT_PIOLIBRARY)
    TARGET_COMPILE_DEFINITIONS(${libName} PUBLIC HOMME_WITHOUT_PIOLIBRARY)
  ENDIF()

  target_link_libraries(${libName} PUBLIC csm_share)
  if (NOT HOMME_BUILD_SCORPIO)
    # Needed for netcdf.mod usage in mesh_mod.F90.
    target_link_libraries(${libName} PUBLIC piof)
  endif()

  IF (CXXLIB_SUPPORTED_CACHE)
    MESSAGE(STATUS "   Linking Fortran with -cxxlib")
    TARGET_LINK_LIBRARIES(${libName} PUBLIC -cxxlib)
  ENDIF ()

  STRING(TOUPPER "${PERFORMANCE_PROFILE}" PERF_PROF_UPPER)
  IF ("${PERF_PROF_UPPER}" STREQUAL "VTUNE")
    TARGET_LINK_LIBRARIES(${libName} PUBLIC ittnotify)
  ENDIF ()

  # COMPOSE_LIBRARY is empty if Compose SL transport is not enabled.
  TARGET_LINK_LIBRARIES(${libName} PUBLIC timing ${COMPOSE_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  IF (HOMME_USE_KOKKOS)
    TARGET_LINK_LIBRARIES(${libName} PUBLIC Kokkos::kokkos)
  ENDIF ()

  IF (HOMME_USE_MKL)
    IF (MKL_TYPE STREQUAL "oneMKL")
      TARGET_LINK_LIBRARIES(${libName} PUBLIC -qmkl)
    ELSEIF (MKL_TYPE STREQUAL "Intel MKL")
      TARGET_LINK_LIBRARIES(${libName} PUBLIC -mkl)
    ENDIF ()
  ELSE()
    IF (NOT HOMME_FIND_BLASLAPACK)
      TARGET_LINK_LIBRARIES(${libName} PUBLIC lapack blas)
    ENDIF()
  ENDIF()

  IF (HAVE_EXTRAE)
    TARGET_LINK_LIBRARIES(${libName} PUBLIC ${Extrae_LIBRARY})
  ENDIF ()

  IF (HOMME_USE_TRILINOS)
    TARGET_LINK_LIBRARIES(${libName} PUBLIC ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  ENDIF()

  IF (HOMME_USE_ARKODE AND "${execType}" STREQUAL "theta-l")
    TARGET_LINK_LIBRARIES(${libName} PUBLIC sundials_farkode)
    TARGET_LINK_LIBRARIES(${libName} PUBLIC sundials_arkode)
    TARGET_LINK_LIBRARIES(${libName} PUBLIC sundials_nvecserial)
    TARGET_LINK_LIBRARIES(${libName} PUBLIC sundials_fnvecserial)
  ENDIF ()

endmacro(createExecLib)

macro(testQuadPrec HOMME_QUAD_PREC)

  TRY_COMPILE(COMPILE_RESULT_VAR
              ${HOMME_BINARY_DIR}/tests/compilerTests/
              ${HOMME_SOURCE_DIR}/cmake/compilerTests/quadTest.f90
              OUTPUT_VARIABLE COMPILE_OUTPUT)

  IF (${COMPILE_RESULT_VAR})
    SET (HOMME_QUAD_PREC TRUE)
    MESSAGE(STATUS "Quadruple-precision supported enabling")
  ELSE ()
    SET (HOMME_QUAD_PREC FALSE)
    MESSAGE(STATUS "Quadruple-precision not supported")
  ENDIF ()
endmacro(testQuadPrec)


macro(setCustomCompilerFlags CUSTOM_FLAGS_FILE SRCS_ALL)

  # Locally reset the compiler flags
  #   This only changes the flags for preqx
  #   these variables get reset outside of this subdir
  SET(CMAKE_Fortran_FLAGS_ORIG "${CMAKE_Fortran_FLAGS}")
  SET(CMAKE_Fortran_FLAGS "")

  # Need to put the OpenMP flag back on the compile line since it
  #   is used for linking
  IF (${ENABLE_OPENMP})
    SET(CMAKE_Fortran_FLAGS "${OpenMP_Fortran_FLAGS}")
  ENDIF ()

  # Need to put the -mmic flag back on the compile line since it
  #   is used for linking
  IF (${ENABLE_INTEL_PHI})
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${INTEL_PHI_FLAGS}")
  ENDIF ()

  # This file should declare the list of files to be exclude
  #   from default compilation and declare compiler options for them
  INCLUDE(${${CUSTOM_FLAGS_FILE}})

  # Remove the custom files from the list of all files
  FOREACH (CUSTOM_FILE ${CUSTOM_FLAG_FILES})
    MESSAGE(STATUS "Applying custom compiler flags to ${CUSTOM_FILE}")
    GET_SOURCE_FILE_PROPERTY(THIS_CUSTOM_FLAGS ${CUSTOM_FILE} COMPILE_FLAGS)
    MESSAGE(STATUS "  ${THIS_CUSTOM_FLAGS}")
    LIST(REMOVE_ITEM ${SRCS_ALL} ${CUSTOM_FILE})
  ENDFOREACH()

  # Compile the rest of the files with the original flags
  SET_SOURCE_FILES_PROPERTIES(${${SRCS_ALL}} PROPERTIES COMPILE_FLAGS
                              "${CMAKE_Fortran_FLAGS_ORIG}")

  # Add the custom files back in to the list of all files
  SET(${SRCS_ALL} ${${SRCS_ALL}} ${CUSTOM_FLAG_FILES})

endmacro(setCustomCompilerFlags)

macro(setSrcMods SRCMODS_PATH SRCS_ALL)
  #Use the files declared in a single SOURCEMODS path as replacements.
  FILE(GLOB SRCMODS_LIST "${${SRCMODS_PATH}}/*.F90")
  SET(TMP_SRCS_ALL ${${SRCS_ALL}})
  FOREACH (SRCMODS_FILE ${SRCMODS_LIST})
   GET_FILENAME_COMPONENT(SRCMOD_NAME ${SRCMODS_FILE} NAME)
   FOREACH (SRCS_ALL_FILE ${TMP_SRCS_ALL})
     GET_FILENAME_COMPONENT(SRCS_ALL_NAME ${SRCS_ALL_FILE} NAME)
     IF (${SRCMOD_NAME} STREQUAL ${SRCS_ALL_NAME})
       LIST(REMOVE_ITEM TMP_SRCS_ALL ${SRCS_ALL_FILE})
       LIST(APPEND TMP_SRCS_ALL ${SRCMODS_FILE})
     ENDIF()
   ENDFOREACH()
  ENDFOREACH()
  SET(${SRCS_ALL} ${TMP_SRCS_ALL} CACHE INTERNAL "modified sources")
endmacro(setSrcMods)
