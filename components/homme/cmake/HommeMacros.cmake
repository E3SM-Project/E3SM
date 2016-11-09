# String representation of pound since it is a cmake comment char.
STRING(ASCII 35 POUND)

function (prc var)
  message("${var} ${${var}}")
endfunction ()

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
  SET(THIS_CONFIG_HC ${CMAKE_CURRENT_BINARY_DIR}/config.h.c)
  SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)

  # First configure the file (which formats the file as C)
  CONFIGURE_FILE(${HOMME_SOURCE_DIR}/src/${execType}/config.h.cmake.in ${THIS_CONFIG_HC})

  # Next reformat the file as Fortran by appending comment lines with an exclamation mark
  EXECUTE_PROCESS(COMMAND sed "s;^/;!/;g"
                     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                     INPUT_FILE ${THIS_CONFIG_HC}
                     OUTPUT_FILE ${THIS_CONFIG_H})

  ADD_DEFINITIONS(-DHAVE_CONFIG_H)
  
  ADD_EXECUTABLE(${execName} ${EXEC_SOURCES})

  # Add this executable to a list 
  SET(EXEC_LIST ${EXEC_LIST} ${execName} CACHE INTERNAL "List of configured executables")

  TARGET_LINK_LIBRARIES(${execName} pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  # Move the module files out of the way so the parallel build 
  # doesn't have a race condition
  SET_TARGET_PROPERTIES(${execName} 
                        PROPERTIES Fortran_MODULE_DIRECTORY ${EXEC_MODULE_DIR})

  IF (HOMME_USE_MKL)
    TARGET_LINK_LIBRARIES(${execName})
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

  INSTALL(TARGETS ${execName} RUNTIME DESTINATION tests)

endmacro(createTestExec)



macro (copyDirFiles testDir)
  # Copy all of the files into the binary dir
  FOREACH (singleFile ${NAMELIST_FILES}) 
    FILE(COPY ${singleFile} DESTINATION ${testDir})
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
  FOREACH (singleFile ${REFSOLN_FILES}) 
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
  SET(THIS_BASELINE_TEST_DIR ${CMAKE_BINARY_DIR}/tests/baseline/${TEST_NAME})

  copyDirFiles(${TEST_DIR})
  copyDirFiles(${THIS_BASELINE_TEST_DIR})

  # Create a run script
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
  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
  SET (TEST_INDEX 1)
  FOREACH (singleFile ${NAMELIST_FILES}) 
    FILE(APPEND ${THIS_TEST_SCRIPT} "TEST_${TEST_INDEX}=\"${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}/${EXEC_NAME} < ${singleFile}\"\n")
    FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
    MATH(EXPR TEST_INDEX "${TEST_INDEX} + 1")
  ENDFOREACH ()
  MATH(EXPR TEST_INDEX "${TEST_INDEX} - 1")
  FILE(APPEND ${THIS_TEST_SCRIPT} "NUM_TESTS=${TEST_INDEX}\n") # new line
  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line

  # openMP runs
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
  SET(REFSOLN_FILES)
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
  IF (NOT "${REFSOLN_FILES}" STREQUAL "")
    MESSAGE(STATUS "  refsoln_files=")
    FOREACH (singleFile ${REFSOLN_FILES}) 
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

# Macro to create the individual tests
macro(createTest testFile)

  SET (THIS_TEST_INPUT ${HOMME_SOURCE_DIR}/test/reg_test/run_tests/${testFile})

  resetTestVariables()

  SET(HOMME_ROOT ${HOMME_SOURCE_DIR})

  INCLUDE(${THIS_TEST_INPUT})

  FILE(GLOB NAMELIST_FILES ${NAMELIST_FILES})
  FILE(GLOB VCOORD_FILES ${VCOORD_FILES})
  FILE(GLOB NCL_FILES ${NCL_FILES})
  FILE(GLOB REFSOLN_FILES ${REFSOLN_FILES})
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
    SET(THIS_TEST_DIR ${CMAKE_BINARY_DIR}/tests/${TEST_NAME})

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
               COMMAND ${CMAKE_BINARY_DIR}/tests/diff_output.sh ${TEST_NAME})

      SET_TESTS_PROPERTIES(${THIS_TEST} PROPERTIES DEPENDS submitAndRunTests)


    ELSE ()

      SET(THIS_TEST "${TEST_NAME}")

      # When not run through a queue each run is ran and then diffed. This is handled by 
      #  the submit_tests.sh script 
      ADD_TEST(NAME ${THIS_TEST} 
               COMMAND ${CMAKE_BINARY_DIR}/tests/submit_tests.sh "${THIS_TEST_RUN_SCRIPT}" "${TEST_NAME}"
               DEPENDS ${EXEC_NAME})

    ENDIF ()

    # Force cprnc to be built when the individual test is run
    SET_TESTS_PROPERTIES(${THIS_TEST} PROPERTIES DEPENDS cprnc)

    # Individual target to rerun and diff the tests
    SET(THIS_TEST_INDIV "test-${TEST_NAME}")

    ADD_CUSTOM_TARGET(${THIS_TEST_INDIV}
             COMMAND ${CMAKE_BINARY_DIR}/tests/submit_tests.sh "${THIS_TEST_RUN_SCRIPT}" "${TEST_NAME}")

    ADD_DEPENDENCIES(${THIS_TEST_INDIV} ${EXEC_NAME})

    # Check target 
    ADD_DEPENDENCIES(check ${EXEC_NAME})

    # Baseline target
    ADD_DEPENDENCIES(baseline ${EXEC_NAME})

    # Force cprnc to be built when the individual test is run
    ADD_DEPENDENCIES(${THIS_TEST_INDIV} cprnc)

    # Now make the Individual targets
    #ADD_CUSTOM_COMMAND(TARGET ${THIS_TEST_INDIV}
    #                   COMMENT "Running the HOMME regression test: ${THIS_TEST}"
    #                   POST_BUILD COMMAND ${CMAKE_CTEST_COMMAND} ARGS --output-on-failure -R ${THIS_TEST_INDIV} 
    #                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  ENDIF ()
endmacro(createTest)

macro(createTests testList)
 
  FOREACH (test ${${testList}})
    createTest(${test})
  ENDFOREACH ()
endmacro(createTests)


macro(testQuadPrec HOMME_QUAD_PREC)

  TRY_COMPILE(COMPILE_RESULT_VAR
              ${CMAKE_BINARY_DIR}/tests/compilerTests/
              ${CMAKE_SOURCE_DIR}/cmake/compilerTests/quadTest.f90
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
