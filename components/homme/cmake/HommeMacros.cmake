# String representation of pound since it is a cmake comment char.
STRING(ASCII 35 POUND)

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

    IF (${CONFIG_FILE_C_STR} STREQUAL ${CONFIG_FILE_C_TMP_STR})
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
  SET_TARGET_PROPERTIES(${execName} PROPERTIES LINKER_LANGUAGE Fortran)
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
    link_to_kokkos(${execName})
  ENDIF ()

  # Move the module files out of the way so the parallel build
  # doesn't have a race condition
  SET_TARGET_PROPERTIES(${execName}
                        PROPERTIES Fortran_MODULE_DIRECTORY ${EXEC_MODULE_DIR})

  IF (HOMME_USE_MKL)
    TARGET_COMPILE_OPTIONS(${execName} PUBLIC -mkl)
    TARGET_LINK_LIBRARIES(${execName} -mkl)
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

  target_link_libraries(${execName} csm_share)
  if (NOT HOMME_BUILD_SCORPIO)
    # Needed for netcdf.mod usage in mesh_mod.F90.
    target_link_libraries(${execName} piof)
  endif()

  IF (CXXLIB_SUPPORTED_CACHE)
    MESSAGE(STATUS "   Linking Fortran with -cxxlib")
    TARGET_LINK_LIBRARIES(${libName} -cxxlib)
  ENDIF ()

  STRING(TOUPPER "${PERFORMANCE_PROFILE}" PERF_PROF_UPPER)
  IF ("${PERF_PROF_UPPER}" STREQUAL "VTUNE")
    TARGET_LINK_LIBRARIES(${libName} ittnotify)
  ENDIF ()

  # COMPOSE_LIBRARY is empty if Compose SL transport is not enabled.
  TARGET_LINK_LIBRARIES(${libName} timing ${COMPOSE_LIBRARY} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  IF (HOMME_USE_KOKKOS)
    TARGET_LINK_LIBRARIES(${libName} kokkos)
  ENDIF ()

  IF (HOMME_USE_MKL)
    TARGET_COMPILE_OPTIONS(${libName} PUBLIC -mkl)
    TARGET_LINK_LIBRARIES(${libName} -mkl)
  ELSE()
    IF (NOT HOMME_FIND_BLASLAPACK)
      TARGET_LINK_LIBRARIES(${libName} lapack blas)
    ENDIF()
  ENDIF()

  IF (HAVE_EXTRAE)
    TARGET_LINK_LIBRARIES(${libName} ${Extrae_LIBRARY})
  ENDIF ()

  IF (HOMME_USE_TRILINOS)
    TARGET_LINK_LIBRARIES(${libName} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  ENDIF()

  IF (HOMME_USE_ARKODE AND "${execType}" STREQUAL "theta-l")
    TARGET_LINK_LIBRARIES(${libName} sundials_farkode)
    TARGET_LINK_LIBRARIES(${libName} sundials_arkode)
    TARGET_LINK_LIBRARIES(${libName} sundials_nvecserial)
    TARGET_LINK_LIBRARIES(${libName} sundials_fnvecserial)
  ENDIF ()

endmacro(createExecLib)

macro (copyDirFiles testDir)
  # Copy all of the files into the binary dir
  FOREACH (singleFile ${NAMELIST_FILES})
    # Some namelist contain cmake variable, to generate
    # multiple testcases with different ne or ndays,
    # so use CONFIGURE_FILE, to replace variables with values
    GET_FILENAME_COMPONENT(fileName ${singleFile} NAME)
    CONFIGURE_FILE(${singleFile} ${testDir}/${fileName})
  ENDFOREACH ()
  FOREACH (singleFile ${VCOORD_FILES})
    FILE(COPY ${singleFile} DESTINATION ${testDir}/vcoord)
  ENDFOREACH ()
  FOREACH (singleFile ${MESH_FILES})
    FILE(COPY ${singleFile} DESTINATION ${testDir})
  ENDFOREACH ()
  FOREACH (singleFile ${NCL_FILES})
    FILE(COPY ${singleFile} DESTINATION ${testDir})
  ENDFOREACH ()
  FOREACH (singleFile ${MESH_FILES})
    FILE(COPY ${singleFile} DESTINATION ${testDir})
  ENDFOREACH ()
  FOREACH (singleFile ${TRILINOS_XML_FILE})
    FILE(COPY ${singleFile} DESTINATION ${testDir})
  ENDFOREACH ()




  # Need to create the movie directory for output
  EXECUTE_PROCESS(COMMAND mkdir -p ${testDir}/movies
    RESULT_VARIABLE Homme_result
    OUTPUT_VARIABLE Homme_output
    ERROR_VARIABLE Homme_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Need to create the restart directory for restart files
  EXECUTE_PROCESS(COMMAND mkdir -p ${testDir}/restart
    RESULT_VARIABLE Homme_result
    OUTPUT_VARIABLE Homme_output
    ERROR_VARIABLE Homme_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

endmacro (copyDirFiles)

macro (setUpTestDir TEST_DIR)

  SET(TEST_DIR_LIST ${TEST_DIR_LIST} ${TEST_DIR})
  SET(THIS_BASELINE_TEST_DIR ${HOMME_BINARY_DIR}/tests/baseline/${TEST_NAME})

  copyDirFiles(${TEST_DIR})
  copyDirFiles(${THIS_BASELINE_TEST_DIR})

  # Create a template TEST.sh for a run script
  SET(THIS_TEST_SCRIPT ${TEST_DIR}/${TEST_NAME}.sh)

  FILE(WRITE  ${THIS_TEST_SCRIPT} "${POUND}!/bin/bash\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND} Test script generated by CMake\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND} Used for running ${TEST_NAME}\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")

  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
  IF ("${NUM_CPUS}" STREQUAL "")
    MESSAGE(FATAL_ERROR "In test ${TEST_NAME} NUM_CPUS not defined. Quitting")
  ENDIF ()
  IF (NOT ${HOMME_QUEUING})
    IF (NOT ${USE_NUM_PROCS} STREQUAL "")
    #IF (USE_NUM_PROCS)
      #FILE(APPEND ${THIS_TEST_SCRIPT} "num_cpus=${USE_NUM_PROCS}\n") # new line
      SET(NUM_CPUS ${USE_NUM_PROCS})
    ELSEIF (${NUM_CPUS} GREATER ${MAX_NUM_PROCS})
      ##MESSAGE(STATUS "For ${TEST_NAME} the requested number of CPU processes is larger than the number available")
      ##MESSAGE(STATUS "  Changing NUM_CPU from ${NUM_CPUS} to ${MAX_NUM_PROCS}")
      ##SET(NUM_CPUS ${MAX_NUM_PROCS})
      #FILE(APPEND ${THIS_TEST_SCRIPT} "num_cpus=${MAX_NUM_PROCS}\n") # new line
      SET(NUM_CPUS ${MAX_NUM_PROCS})
    ENDIF ()
  ENDIF ()
  FILE(APPEND ${THIS_TEST_SCRIPT} "NUM_CPUS=${NUM_CPUS}\n") # new line
  FILE(APPEND ${THIS_TEST_SCRIPT} "EXEC=${EXEC_NAME}\n")

  #not making it more general due to jsrun, p9/gpu, loading modules...
  if (HOMME_MACHINE MATCHES "summit-gpu")
    #cmake 3.17
    #foreach(varn varv IN ZIP_LISTS ${varnames} ${varvals})
    #  file(APPEND ${THIS_TEST_SCRIPT} "${varn}=\"${varv}\"")
    #  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
    #endforeach()

    list(LENGTH varnames len1)
    math(EXPR len2 "${len1} - 1")

    foreach(val RANGE ${len2})
      list(GET varnames ${val} varn)
      list(GET varvals ${val} varv)
      file(APPEND ${THIS_TEST_SCRIPT} "${varn}=\"${varv}\"")
      FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
    endforeach()
  endif()

  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line

  SET (TEST_INDEX 1)
  FOREACH (singleFile ${NAMELIST_FILES})
    GET_FILENAME_COMPONENT(fileName ${singleFile} NAME)
    FILE(APPEND ${THIS_TEST_SCRIPT} "TEST_${TEST_INDEX}=\"${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}/${EXEC_NAME} < ${fileName}\"\n")
    FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
    MATH(EXPR TEST_INDEX "${TEST_INDEX} + 1")
  ENDFOREACH ()
  MATH(EXPR TEST_INDEX "${TEST_INDEX} - 1")
  FILE(APPEND ${THIS_TEST_SCRIPT} "NUM_TESTS=${TEST_INDEX}\n") # new line
  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line

  # openMP runs

  IF (EXEC_NAME MATCHES "kokkos$")
    #not doing logic around omp_num_mpi...
    FILE(APPEND ${THIS_TEST_SCRIPT} "OMP_NUMBER_THREADS_KOKKOS=${OMP_NUM_THREADS}\n") # new line
  ENDIF()

  IF (NOT "${OMP_NAMELIST_FILES}" STREQUAL "")
    IF (${ENABLE_HORIZ_OPENMP})
      FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
      FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND} OpenMP Tests\n")
      FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
      FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
      MATH(EXPR OMP_MOD "${NUM_CPUS} % ${OMP_NUM_THREADS}")
      IF (NOT ${OMP_MOD} EQUAL 0)
      # MT: why is this a fatal error?  i'm just trying to build preqx and not run regression tests.
        MESSAGE(STATUS "In test ${TEST_NAME} NUM_CPUS not divisible by OMP_NUM_THREADS. Quitting.")
      ENDIF ()
      MATH(EXPR OMP_NUM_MPI "${NUM_CPUS} / ${OMP_NUM_THREADS}")
      FILE(APPEND ${THIS_TEST_SCRIPT} "OMP_NUM_MPI=${OMP_NUM_MPI}\n") # new line
      FILE(APPEND ${THIS_TEST_SCRIPT} "OMP_NUMBER_THREADS=${OMP_NUM_THREADS}\n") # new line
      FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
      SET (TEST_INDEX 1)
      FOREACH (singleFile ${OMP_NAMELIST_FILES})
        FILE(APPEND ${THIS_TEST_SCRIPT} "OMP_TEST_${TEST_INDEX}=\"${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}/${EXEC_NAME} < ${singleFile}\"\n")
        FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
        MATH(EXPR TEST_INDEX "${TEST_INDEX} + 1")
      ENDFOREACH ()
      MATH(EXPR TEST_INDEX "${TEST_INDEX} - 1")
      FILE(APPEND ${THIS_TEST_SCRIPT} "OMP_NUM_TESTS=${TEST_INDEX}\n") # new line
    ELSE ()
      MESSAGE(STATUS "  Not including OpenMP tests")
    ENDIF()
  ENDIF ()

  # Add this test to the list of tests
  MATH(EXPR NUM_TEST_FILES "${NUM_TEST_FILES} + 1")
  FILE (APPEND ${HOMME_TEST_LIST} "TEST_FILE_${NUM_TEST_FILES}=${THIS_TEST_SCRIPT}\n")

  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line

  # Deal with the Netcdf output files
  FILE(APPEND ${THIS_TEST_SCRIPT} "NC_OUTPUT_FILES=\"")
  FOREACH (singleFile ${NC_OUTPUT_FILES})
    FILE(APPEND ${THIS_TEST_SCRIPT} "${singleFile} ")
  ENDFOREACH ()
  # Add the OPENMP netcdf outpuf files
  IF (${ENABLE_HORIZ_OPENMP})
    FOREACH (singleFile ${OMP_NC_OUTPUT_FILES})
      FILE(APPEND ${THIS_TEST_SCRIPT} "${singleFile} ")
    ENDFOREACH ()
  ENDIF()
  FILE(APPEND ${THIS_TEST_SCRIPT} "\"\n")

  FILE(APPEND ${THIS_TEST_SCRIPT} "NC_OUTPUT_REF=\"")
  FOREACH (singleFile ${NC_OUTPUT_REF})
    FILE(APPEND ${THIS_TEST_SCRIPT} "${singleFile} ")
  ENDFOREACH ()
  FILE(APPEND ${THIS_TEST_SCRIPT} "\"\n")

  FILE(APPEND ${THIS_TEST_SCRIPT} "NC_OUTPUT_CHECKREF=\"")
  FOREACH (singleFile ${NC_OUTPUT_CHECKREF})
    FILE(APPEND ${THIS_TEST_SCRIPT} "${singleFile} ")
  ENDFOREACH ()
  FILE(APPEND ${THIS_TEST_SCRIPT} "\"\n")

  FILE(APPEND ${THIS_TEST_SCRIPT} "TESTCASE_REF_TOL=\"")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${TESTCASE_REF_TOL}")
  FILE(APPEND ${THIS_TEST_SCRIPT} "\"\n")



endmacro (setUpTestDir)

macro(resetTestVariables)
  # Reset the variables
  SET(TEST_NAME)
  SET(EXEC_NAME)
  SET(NAMELIST_FILES)
  SET(VCOORD_FILES)
  SET(MESH_FILES)
  SET(NCL_FILES)
  SET(MESH_FILES)
  SET(NC_OUTPUT_FILES)
  SET(OMP_NC_OUTPUT_FILES)
  SET(NC_OUTPUT_REF)
  SET(NC_OUTPUT_CHECKREF)
  SET(TESTCASE_REF_TOL)
  SET(NUM_CPUS)
  SET(OMP_NAMELIST_FILES)
  SET(OMP_NUM_THREADS)
  SET(TRILINOS_XML_FILE)
endmacro(resetTestVariables)

macro(printTestSummary)
  MESSAGE(STATUS "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  MESSAGE(STATUS "Summary of test ${TEST_NAME}")
  IF (NOT "${NAMELIST_FILES}" STREQUAL "")
    MESSAGE(STATUS "  namelist_files=")
    FOREACH (singleFile ${NAMELIST_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${VCOORD_FILES}" STREQUAL "")
    MESSAGE(STATUS "  vcoord_files=")
    FOREACH (singleFile ${VCOORD_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${MESH_FILES}" STREQUAL "")
    MESSAGE(STATUS "  mesh_files=")
    FOREACH (singleFile ${MESH_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${NCL_FILES}" STREQUAL "")
    MESSAGE(STATUS "  ncl_files=")
    FOREACH (singleFile ${NCL_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${NC_OUTPUT_FILES}" STREQUAL "")
    MESSAGE(STATUS "  nc_output_files=")
    FOREACH (singleFile ${NC_OUTPUT_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${NC_OUTPUT_REF}" STREQUAL "")
    MESSAGE(STATUS "  nc_output_files=")
    FOREACH (singleFile ${NC_OUTPUT_REF})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${NC_OUTPUT_CHECKREF}" STREQUAL "")
    MESSAGE(STATUS "  nc_output_files=")
    FOREACH (singleFile ${NC_OUTPUT_CHECKREF})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${TESTCASE_REF_TOL}" STREQUAL "")
    MESSAGE(STATUS "  testcase_ref_tol=")
    MESSAGE(STATUS "  ${TESTCASE_REF_TOL}")
  ENDIF ()
  IF (NOT "${MESH_FILES}" STREQUAL "")
    MESSAGE(STATUS "  mesh_files=")
    FOREACH (singleFile ${MESH_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${OMP_NAMELIST_FILES}" STREQUAL "")
    MESSAGE(STATUS "  omp_namelist_files=")
    FOREACH (singleFile ${OMP_NAMELIST_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${OMP_NC_OUTPUT_FILES}" STREQUAL "")
    MESSAGE(STATUS "  omp_nc_output_files=")
    FOREACH (singleFile ${OMP_NC_OUTPUT_FILES})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()
  IF (NOT "${TRILINOS_XML_FILE}" STREQUAL "")
    MESSAGE(STATUS "  trilinos_xml_file=")
    FOREACH (singleFile ${TRILINOS_XML_FILE})
      MESSAGE(STATUS "    ${singleFile}")
    ENDFOREACH ()
  ENDIF ()

endmacro(printTestSummary)

macro (set_homme_tests_parameters testFile profile)
#example on how to customize tests lengths
#  if ("${testFile}" MATCHES ".*-moist.*")

#theta tests for kokkos targets
#HV scaled at 3.0 rate: 
#ne30 1e15, ne2 3.4e18, ne4 4.2e17
#ne6 1.25e17, ne13 1.23e16, ne1024 2.5e10
  if ("${testFile}" MATCHES "theta-f")
    if ("${profile}" STREQUAL "dev")
      set (HOMME_TEST_NE 2)
      set (HOMME_TEST_NDAYS 1)
      set (HOMME_TEST_NU 3.4e18)
    elseif ("${profile}" STREQUAL "short")
      set (HOMME_TEST_NE 6)
      set (HOMME_TEST_NDAYS 1) #should be 9 but for GB reduce
      set (HOMME_TEST_NU 1.25e17)
    else ()
      set (HOMME_TEST_NE 13)
      set (HOMME_TEST_NDAYS 1) #should be 9 but for GB reduce
      set (HOMME_TEST_NU 1.2e16)
    endif ()
#preqx tests for kokkos targets
  elseif ("${testFile}" MATCHES "preqx-nlev")
    if ("${testFile}" MATCHES ".*-tensorhv-*")
      if ("${profile}" STREQUAL "dev")
        set (HOMME_TEST_NE 2)
        set (HOMME_TEST_NDAYS 1)
      elseif ("${profile}" STREQUAL "short")
        set (HOMME_TEST_NE 4)
        set (HOMME_TEST_NDAYS 3)
      else ()
        set (HOMME_TEST_NE 12)
        set (HOMME_TEST_NDAYS 3)
      endif ()
    else ()
      if ("${profile}" STREQUAL "dev")
        set (HOMME_TEST_NE 2)
        set (HOMME_TEST_NDAYS 1)
      elseif ("${profile}" STREQUAL "short")
        set (HOMME_TEST_NE 4)
        set (HOMME_TEST_NDAYS 9)
      else ()
        set (HOMME_TEST_NE 12)
        set (HOMME_TEST_NDAYS 9)
      endif ()
    endif()
  endif ()

endmacro ()

# Macro to create the individual tests
macro(createTest testFile)

  SET (THIS_TEST_INPUT ${HOMME_SOURCE_DIR}/test/reg_test/run_tests/${testFile})

  resetTestVariables()

  SET(HOMME_ROOT ${HOMME_SOURCE_DIR})

  INCLUDE(${THIS_TEST_INPUT})

  IF (DEFINED PROFILE)
    set_homme_tests_parameters(${testFile} ${PROFILE})
    if ("${TEST_NAME}" MATCHES "theta-f")
      set (TEST_NAME "${TEST_NAME}-ne${HOMME_TEST_NE}-nu${HOMME_TEST_NU}-ndays${HOMME_TEST_NDAYS}")
    elseif ("${TEST_NAME}" MATCHES "preqx-nlev")
      set (TEST_NAME "${TEST_NAME}-ne${HOMME_TEST_NE}-ndays${HOMME_TEST_NDAYS}")
    endif()
  ENDIF ()

  FILE(GLOB NAMELIST_FILES ${NAMELIST_FILES})
  FILE(GLOB VCOORD_FILES ${VCOORD_FILES})
  FILE(GLOB NCL_FILES ${NCL_FILES})
  FILE(GLOB MESH_FILES ${MESH_FILES})
  FILE(GLOB OMP_NAMELIST_FILES ${OMP_NAMELIST_FILES})
  FILE(GLOB TRILINOS_XML_FILE ${TRILINOS_XML_FILE})

  # Determine if the executable this tests depeds upon is built
  LIST(FIND EXEC_LIST ${EXEC_NAME} FIND_INDEX)

  IF (${FIND_INDEX} LESS 0)
    MESSAGE(STATUS "Not configuring test ${TEST_NAME} since it depends upon the executable ${EXEC_NAME}
                    which isn't built with this configuration")
  ELSE ()

    MESSAGE(STATUS "Adding test: ${TEST_NAME}, using exec ${EXEC_NAME}")

    OPTION(TEST_SUMMARY "Print out information about the tests" OFF)

    IF (${TEST_SUMMARY})
      printTestSummary()
    ENDIF()

    # Set up the directory
    SET(THIS_TEST_DIR ${HOMME_BINARY_DIR}/tests/${TEST_NAME})

    # Set up the test directory for both the baseline and the comparison tests
    setUpTestDir(${THIS_TEST_DIR})
    #setUpTestDir(${THIS_BASELINE_TEST_DIR})

    # The test (not the baseline) run script
    SET(THIS_TEST_RUN_SCRIPT "${THIS_TEST_DIR}/${TEST_NAME}-run.sh")

    IF (${HOMME_QUEUING})

      SET(THIS_TEST "${TEST_NAME}-diff")

      # When run through the queue the runs are submitted and ran in
      #   submitAndRunTests, and diffed in the subsequent tests
      ADD_TEST(NAME ${THIS_TEST}
               COMMAND ${HOMME_BINARY_DIR}/tests/diff_output.sh ${TEST_NAME})

      SET_TESTS_PROPERTIES(${THIS_TEST} PROPERTIES DEPENDS submitAndRunTests)


    ELSE ()

      SET(THIS_TEST "${TEST_NAME}")

      # When not run through a queue each run is ran and then diffed. This is handled by
      #  the submit_tests.sh script
      ADD_TEST(NAME ${THIS_TEST}
               COMMAND ${HOMME_BINARY_DIR}/tests/submit_tests.sh "${THIS_TEST_RUN_SCRIPT}" "${TEST_NAME}"
               DEPENDS ${EXEC_NAME})

    ENDIF ()

    # Force cprnc to be built when the individual test is run
    SET_TESTS_PROPERTIES(${THIS_TEST} PROPERTIES DEPENDS cprnc)

    IF (DEFINED PROFILE)
      SET_TESTS_PROPERTIES(${THIS_TEST} PROPERTIES LABELS ${PROFILE})
    ENDIF()

    # Individual target to rerun and diff the tests
    SET(THIS_TEST_INDIV "test-${TEST_NAME}")

    ADD_CUSTOM_TARGET(${THIS_TEST_INDIV}
             COMMAND ${HOMME_BINARY_DIR}/tests/submit_tests.sh "${THIS_TEST_RUN_SCRIPT}" "${TEST_NAME}")

    ADD_DEPENDENCIES(${THIS_TEST_INDIV} ${EXEC_NAME})

    # test execs
    ADD_DEPENDENCIES(test-execs ${EXEC_NAME})

    # Check target
    ADD_DEPENDENCIES(check ${EXEC_NAME})

    # Baseline target
    ADD_DEPENDENCIES(baseline ${EXEC_NAME})

    # Force cprnc to be built when the individual test is run
    IF (TARGET cprnc)
      ADD_DEPENDENCIES(${THIS_TEST_INDIV} cprnc)
    ENDIF()

    # This helped in some builds on GPU, where the test hanged for a VERY long time
    IF (NOT "${TIMEOUT}" STREQUAL "")
      SET_TESTS_PROPERTIES(${THIS_TEST} PROPERTIES TIMEOUT ${TIMEOUT})
    ENDIF ()

    # Now make the Individual targets
    #ADD_CUSTOM_COMMAND(TARGET ${THIS_TEST_INDIV}
    #                   COMMENT "Running the HOMME regression test: ${THIS_TEST}"
    #                   POST_BUILD COMMAND ${CMAKE_CTEST_COMMAND} ARGS --output-on-failure -R ${THIS_TEST_INDIV}
    #                   WORKING_DIRECTORY ${HOMME_BINARY_DIR})

  ENDIF ()
endmacro(createTest)

macro (createTests testList)
  FOREACH (test ${${testList}})
    createTest(${test})
  ENDFOREACH ()
endmacro (createTests)

macro(createTestsWithProfile testList testsProfile)
  SET (PROFILE ${testsProfile})
  FOREACH (test ${${testList}})
    createTest(${test})
  ENDFOREACH ()
  UNSET(PROFILE)
endmacro(createTestsWithProfile)

# Make a list of all testing profiles no more intensive than the given profile.
function (make_profiles_up_to profile profiles)
  string (TOLOWER "${profile}" profile_ci)
  set (tmp)
  if ("${profile_ci}" STREQUAL "dev")
    list (APPEND tmp "dev")
  elseif ("${profile_ci}" STREQUAL "short")
    list (APPEND tmp "dev")
    list (APPEND tmp "short")

  elseif ("${profile_ci}" STREQUAL "nightly")
    list (APPEND tmp "dev")
    list (APPEND tmp "short")
    list (APPEND tmp "nightly")
  else ()
    message (FATAL_ERROR "Testing profile '${profile}' not implemented.")
  endif ()
  set (profiles "${tmp}" PARENT_SCOPE)
endfunction ()

MACRO(CREATE_CXX_VS_F90_TESTS_WITH_PROFILE TESTS_LIST testProfile)

  FOREACH (TEST ${${TESTS_LIST}})
    SET (TEST_FILE_F90 "${TEST}.cmake")

    set_homme_tests_parameters(${TEST} ${testProfile})
    set (PROFILE ${testProfile})
    INCLUDE (${HOMME_SOURCE_DIR}/test/reg_test/run_tests/${TEST_FILE_F90})

    if ("${TEST}" MATCHES "theta-f")
      SET (TEST_NAME_SUFFIX "ne${HOMME_TEST_NE}-nu${HOMME_TEST_NU}-ndays${HOMME_TEST_NDAYS}")
    elseif ("${TEST}" MATCHES "preqx-nlev")
      SET (TEST_NAME_SUFFIX "ne${HOMME_TEST_NE}-ndays${HOMME_TEST_NDAYS}")
    endif ()

    SET (F90_TEST_NAME "${TEST}-${TEST_NAME_SUFFIX}")
    SET (CXX_TEST_NAME "${TEST}-kokkos-${TEST_NAME_SUFFIX}")
    SET (F90_DIR ${HOMME_BINARY_DIR}/tests/${F90_TEST_NAME})
    SET (CXX_DIR ${HOMME_BINARY_DIR}/tests/${CXX_TEST_NAME})

    # Compare netcdf output files bit-for-bit AND compare diagnostic lines
    # in the raw output files
    SET (TEST_NAME "${TEST}-${TEST_NAME_SUFFIX}_cxx_vs_f90")
    MESSAGE ("-- Creating cxx-f90 comparison test ${TEST_NAME}")

    CONFIGURE_FILE (${HOMME_SOURCE_DIR}/cmake/CxxVsF90.cmake.in
                    ${HOMME_BINARY_DIR}/tests/${CXX_TEST_NAME}/CxxVsF90.cmake @ONLY)

    ADD_TEST (NAME "${TEST_NAME}"
              COMMAND ${CMAKE_COMMAND} -P CxxVsF90.cmake
              WORKING_DIRECTORY ${HOMME_BINARY_DIR}/tests/${CXX_TEST_NAME})

    SET_TESTS_PROPERTIES(${TEST_NAME} PROPERTIES DEPENDS "${F90_TEST_NAME};${CXX_TEST_NAME}"
      LABELS ${testProfile})
  ENDFOREACH ()
ENDMACRO(CREATE_CXX_VS_F90_TESTS_WITH_PROFILE)

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
