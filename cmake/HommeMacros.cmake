
# String representation of pound since it is a cmake comment char.
STRING(ASCII 35 POUND)

# Macro to create the individual tests
macro(createTestExec execName execType execSources macroNP macroPLEV macroUSE_PIO
                     macroWITH_ENERGY)

  # Backup the cmake variables
  SET(tempNP ${NUM_POINTS})
  SET(tempPLEV ${NUM_PLEV})
  SET(tempUSE_PIO ${PIO})
  SET(tempWITH_ENERGY ${ENERGY_DIAGNOSTICS})

  # Set the variable to the macro variables
  SET(NUM_POINTS ${macroNP})
  SET(NUM_PLEV ${macroPLEV})

  IF (${macroUSE_PIO})
    SET(PIO TRUE)
    SET(PIO_INTERP)
  ELSE ()
    SET(PIO)
    SET(PIO_INTERP TRUE)
  ENDIF ()

  SET(ENERGY_DIAGNOSTICS ${macroWITH_ENERGY})

  # This is needed to compile the test executables with the correct options
  SET(THIS_CONFIG_H ${CMAKE_CURRENT_BINARY_DIR}/config.h)
  CONFIGURE_FILE(${HOMME_SOURCE_DIR}/src/${execType}/config.h.cmake.in ${THIS_CONFIG_H})

  ADD_DEFINITIONS(-DHAVE_CONFIG_H)

  ADD_EXECUTABLE(${execName} ${${execSources}})

  # More thought needs to go into these options
  #IF (macroWITH_ENERGY)
  #ENDIF ()

  TARGET_LINK_LIBRARIES(${execName} pio timing ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

  # Move the module files out of the way so the parallel build 
  # doesn't have a race condition
  SET_TARGET_PROPERTIES(${execName} 
                        PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${execName}_modules")

  IF (NOT HOMME_FIND_BLASLAPACK)
    TARGET_LINK_LIBRARIES(${execName} lapack blas)
    ADD_DEPENDENCIES(${execName} blas lapack)
  ENDIF()

  INSTALL(TARGETS ${execName} RUNTIME DESTINATION tests)

  # Restore the original the cmake variables
  SET(NUM_POINTS ${tempNP})
  SET(NUM_PLEV ${tempPLEV})
  SET(PIO ${tempUSE_PIO})
  SET(ENERGY_DIAGNOSTICS ${tempWITH_ENERGY})

endmacro(createTestExec)

macro(getValue thisLine thisVar)

  # Turn the line into a list
  STRING(REPLACE "=" ";" LINE_LIST ${thisLine})

  # Get the 0th and 1th items of the list
  #LIST(GET LINE_LIST 0 thisKey)
  LIST(GET LINE_LIST 1 thisVal)

  # Remove whitespace  and set to the variable defined by thisVar
  STRING(STRIP ${thisVal} ${thisVar})
  STRING(REPLACE " " ";" ${thisVar} ${${thisVar}})

endmacro(getValue)

macro(parseTestFile fileName testName execName namelistFiles vcoordFiles
                    nclFiles refsolnFiles numCpus ompNamelistFiles ompNumThreads)

  # Reset the variables
  SET(${testName})
  SET(${execName})
  SET(${namelistFiles})
  SET(${vcoordFiles})
  SET(${nclFiles})
  SET(${refsolnFiles})

  # Need to do something with this
  SET(namelist_dir namelists/little_endian)
  SET(HOMME_ROOT ${HOMME_SOURCE_DIR})

  FILE(STRINGS ${fileName}
       Homme_Raw_Paths
       LIMIT_COUNT 100)

  # Go through each 
  FOREACH(LINE ${Homme_Raw_Paths})

    # Needed to propely concatenate text across multiple lines
    STRING(REPLACE "\"" "" LINE ${LINE})

    IF (NOT ${LINE} MATCHES "^${POUND}")
      STRING(REPLACE "\;" "" REP_LINE ${LINE})
      IF (${REP_LINE} MATCHES "^test_name")
        getValue(${REP_LINE} ${testName})
        SET(test_name ${${testName}})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^exec_name")
        getValue(${REP_LINE} ${execName})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^namelist_files")
        getValue(${REP_LINE} ${namelistFiles})
        FILE(GLOB ${namelistFiles} ${${namelistFiles}})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^vcoord_files")
        getValue(${REP_LINE} ${vcoordFiles})
        FILE(GLOB ${vcoordFiles} ${${vcoordFiles}})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^ncl_files")
        getValue(${REP_LINE} ${nclFiles})
        FILE(GLOB ${nclFiles} ${${nclFiles}})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^refsoln_files")
        getValue(${REP_LINE} ${refsolnFiles})
        FILE(GLOB ${refsolnFiles} ${${refsolnFiles}})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^NUM_CPUS")
        getValue(${REP_LINE} ${numCpus})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^OMP_NUM_THREADS")
        getValue(${REP_LINE} ${ompNumThreads})
      ENDIF ()
      IF (${REP_LINE} MATCHES "^omp_namelist_files")
        getValue(${REP_LINE} ${ompNamelistFiles})
        FILE(GLOB ${ompNamelistFiles} ${${ompNamelistFiles}})
      ENDIF ()
    ENDIF()
  ENDFOREACH()

endmacro(parseTestFile)

# Macro to create the individual tests
macro(createTest testName)

  SET (THIS_TEST_INPUT ${HOMME_SOURCE_DIR}/test/reg_test/run_tests/${testName}.in)

  # Weird this isn't getting reset for somre reason
  SET (OMP_NAMELIST_FILES "")

  parseTestFile(${THIS_TEST_INPUT}
                TEST_NAME EXEC_NAME NAMELIST_FILES VCOORD_FILES NCL_FILES
                REFSOLN_FILES NUM_CPUS OMP_NAMELIST_FILES OMP_NUM_THREADS)

  # Copy the file
  #FILE(COPY ${THIS_TEST_INPUT} DESTINATION ${CMAKE_BINARY_DIR}/tests/${TEST_NAME}/)

  MESSAGE(STATUS "Adding test: ${TEST_NAME}, using exec ${EXEC_NAME}")

  IF (TEST_SUMMARY) 
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
    IF (NOT "${OMP_NAMELIST_FILES}" STREQUAL "")
      MESSAGE(STATUS "  omp_namelist_files=")
      FOREACH (singleFile ${OMP_NAMELIST_FILES}) 
        MESSAGE(STATUS "    ${singleFile}")
      ENDFOREACH () 
    ENDIF ()
  ENDIF()

  # Set up the directory
  SET(THIS_TEST_DIR ${CMAKE_BINARY_DIR}/tests/${TEST_NAME})

  # Copy all of the files into the binary dir
  FOREACH (singleFile ${NAMELIST_FILES}) 
    FILE(COPY ${singleFile} DESTINATION ${THIS_TEST_DIR})
  ENDFOREACH () 
  FOREACH (singleFile ${VCOORD_FILES}) 
    FILE(COPY ${singleFile} DESTINATION ${THIS_TEST_DIR}/vcoord)
  ENDFOREACH () 
  FOREACH (singleFile ${NCL_FILES}) 
    FILE(COPY ${singleFile} DESTINATION ${THIS_TEST_DIR})
  ENDFOREACH () 
  FOREACH (singleFile ${REFSOLN_FILES}) 
    FILE(COPY ${singleFile} DESTINATION ${THIS_TEST_DIR})
  ENDFOREACH () 

  # Create a run script
  SET(THIS_TEST_SCRIPT ${THIS_TEST_DIR}/${TEST_NAME}.sh)

  FILE(WRITE  ${THIS_TEST_SCRIPT} "${POUND}!/bin/bash\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND} Test script generated by CMake\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND} Used for running ${TEST_NAME}\n")
  FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")

  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
  IF (NOT ${NUM_CPUS} STREQUAL "")
    FILE(APPEND ${THIS_TEST_SCRIPT} "num_cpus=${NUM_CPUS}\n") # new line
  ELSE ()
    MESSAGE(FATAL_ERROR "In test ${testName} NUM_CPUS not defined. Quitting")
  ENDIF ()
  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
  SET (TEST_INDEX 1)
  FOREACH (singleFile ${NAMELIST_FILES}) 
    FILE(APPEND ${THIS_TEST_SCRIPT} "test${TEST_INDEX}=\"${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}/${EXEC_NAME} < ${singleFile}\"\n")
    FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
    MATH(EXPR TEST_INDEX "${TEST_INDEX} + 1")
  ENDFOREACH ()
  MATH(EXPR TEST_INDEX "${TEST_INDEX} - 1")
  FILE(APPEND ${THIS_TEST_SCRIPT} "num_tests=${TEST_INDEX}\n") # new line
  FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line

  # openMP runs
  IF (NOT "${OMP_NAMELIST_FILES}" STREQUAL "")
    FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
    FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND} OpenMP Tests\n")
    FILE(APPEND ${THIS_TEST_SCRIPT} "${POUND}===============================\n")
    FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
    MATH(EXPR OMP_MOD "${NUM_CPUS} % ${OMP_NUM_THREADS}")
    IF (NOT ${OMP_MOD} EQUAL 0)
      MESSAGE(FATAL_ERROR "In test ${testName} NUM_CPUS not divisible by OMP_NUM_THREADS. Quitting.")
    ENDIF ()
    MATH(EXPR OMP_NUM_MPI "${NUM_CPUS} / ${OMP_NUM_THREADS}")
    FILE(APPEND ${THIS_TEST_SCRIPT} "omp_num_mpi=${OMP_NUM_MPI}\n") # new line
    FILE(APPEND ${THIS_TEST_SCRIPT} "omp_number_threads=${OMP_NUM_THREADS}\n") # new line
    FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
    SET (TEST_INDEX 1)
    FOREACH (singleFile ${OMP_NAMELIST_FILES}) 
      FILE(APPEND ${THIS_TEST_SCRIPT} "omp_test${TEST_INDEX}=\"${CMAKE_CURRENT_BINARY_DIR}/${EXEC_NAME}/${EXEC_NAME} < ${singleFile}\"\n")
      FILE(APPEND ${THIS_TEST_SCRIPT} "\n") # new line
      MATH(EXPR TEST_INDEX "${TEST_INDEX} + 1")
    ENDFOREACH () 
    MATH(EXPR TEST_INDEX "${TEST_INDEX} - 1")
    FILE(APPEND ${THIS_TEST_SCRIPT} "omp_num_tests=${TEST_INDEX}\n") # new line
  ENDIF ()

  # Need to create the movie directory for output
  EXECUTE_PROCESS(COMMAND mkdir -p ${THIS_TEST_DIR}/movies
    RESULT_VARIABLE Homme_result
    OUTPUT_VARIABLE Homme_output
    ERROR_VARIABLE Homme_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  # Add this test to the list of tests
  MATH(EXPR NUM_TEST_FILES "${NUM_TEST_FILES} + 1")
  FILE (APPEND ${HOMME_TEST_LIST} "test_file${NUM_TEST_FILES}=${THIS_TEST_SCRIPT}\n")

  ADD_TEST(${testName} ${THIS_TEST_SCRIPT})

endmacro(createTest)

macro(testQuadPrec HOMME_QUAD_PREC)
  TRY_RUN(RUN_RESULT_VAR COMPILE_RESULT_VAR
          ${CMAKE_BINARY_DIR}/compiler_test/ ${CMAKE_SOURCE_DIR}/cmake/quadTest.f90
          COMPILE_OUTPUT_VARIABLE comp 
          )

  IF (${COMPILE_RESULT_VAR} AND ${RUN_RESULT_VAR} STREQUAL 0)
    SET (${HOMME_QUAD_PREC} TRUE)
    MESSAGE(STATUS "Quadruple-precision supported enabling")
  ELSE ()
    SET (HOMME_QUAD_PREC FALSE)
    MESSAGE(STATUS "Quadruple-precision not supported")
  ENDIF ()
endmacro(testQuadPrec)
