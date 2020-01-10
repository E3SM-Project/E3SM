##### this is needed if we want dereferencing work as now in test-lists.sh
#CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#add a list of tests in test-list.cmake into theta HOMME_THETA_TESTS_WITH_PROFILE_1 list 
#example of a test name:
#theta-f0-tt5-hvs1-hvst0-r3-qz1-nutopon
#the list will be used to generate corresponding ${test}-kokkos.cmake and ${test}.cmake
#run 'source create-theta-tests-MAIN.sh ' to generate these *cmake files
#then config homme


SET(BUILD_HOMME_THETA_KOKKOS TRUE)
#message("${BUILD_HOMME_THETA_KOKKOS}")
include(test-list.cmake)

#message("${HOMME_THETA_TESTS_WITH_PROFILE_1}")

foreach(jj ${HOMME_THETA_TESTS_WITH_PROFILE_1})

  #message(${jj})
  #parse jj
  string(REPLACE "-" ";" jjlist ${jj})
  #message("${jjlist}")

#most values are assumend to be integers <10, only ttype can be bigger
#see below to parse bigger ints

  #form:
  list(GET jjlist 1 aa)
  #now get a value from aa
  string(REGEX MATCH "[0-9]" bb "${aa}")
  message("adv form is ${bb}")
  set(AAADVFORM "${bb}")

  #ttype:
  list(GET jjlist 2 aa)
  #now get a value from aa
  string(REGEX MATCH "[0-9]+" bb "${aa}")
  message("ttype is ${bb}")
  set(AATTYPE "${bb}")
  #sett hy mode based on ttype
  if("${bb}" STREQUAL "5")
    set(AAHYMODE "true")
  elseif("${bb}" STREQUAL "9")
    set(AAHYMODE "false")
  elseif("${bb}" STREQUAL "10")
    set(AAHYMODE "false")
  else()
    message(FATAL_ERROR "ttype should be 5,9, or 10")
  endif() 

  #hvs:
  list(GET jjlist 3 aa)
  #now get a value from aa
  string(REGEX MATCH "[0-9]" bb "${aa}")
  message("hvs is ${bb}")
  set(AAHVS "${bb}")

  #hvst:
  list(GET jjlist 4 aa)
  #now get a value from aa
  string(REGEX MATCH "[0-9]" bb "${aa}")
  message("hvst is ${bb}")
  set(AATOM "${bb}")

  #rsplit:
  list(GET jjlist 5 aa)
  #now get a value from aa
  string(REGEX MATCH "[0-9]" bb "${aa}")
  message("rsplit is ${bb}")
  set(AARSPLIT "${bb}")

  #qsize:
  list(GET jjlist 6 aa)
  #now get a value from aa
  string(REGEX MATCH "[0-9]" bb "${aa}")
  message("qsize is ${bb}")
  set(AAQSIZE "${bb}")

  #nutop:
  list(GET jjlist 7 aa)
  if("${aa}" STREQUAL "nutopon")
    set(AANTOP 2.5e5)
  elseif("${aa}" STREQUAL "nutopoff")
    set(AANTOP 0)
  else()
    message(FATAL_ERROR "only allowed are nutopon/nutopoff")
  endif()

#since we are moving to change nu based on profile IN THETA tests only
#(cause profiles set NE) we would need to change nudiv based on profile
#easiest way to do it is to introduce a nu_div factor,
#factor = nu/nu_div, and then carry math for nudiv when configure_file(namelist)
#is called. 
#instead, we drop nudiv testing in theta because it is a velocity HV op that
#is tested in preqx. 
#  #nudiv:
#  list(GET jjlist 8 aa)
#  if("${aa}" STREQUAL "nudivon")
#    set(AANDIV 1e15)
#  elseif("${aa}" STREQUAL "nudivoff")
#    set(AANDIV 0)
#  else()
#    message(FATAL_ERROR "only allowed are nudivon/nudivoff")
#  endif()

  #two iterations, with kokkos suffix and without
  #cannot make a list with empty string, doing repetition

  set(AAIFKOKKOS "")
  set(AANAME "${jj}${AAIFKOKKOS}")
  configure_file(create-theta-tests-varlist.sh create-theta-temporary.sh)
  execute_process(COMMAND bash create-theta-tests-generate-cmake-file.sh
    RESULT_VARIABLE RESV
    OUTPUT_VARIABLE OUTV
    ERROR_VARIABLE ERRV)

  set(AAIFKOKKOS "-kokkos")
  set(AANAME "${jj}${AAIFKOKKOS}")
  configure_file(create-theta-tests-varlist.sh create-theta-temporary.sh)
  execute_process(COMMAND bash create-theta-tests-generate-cmake-file.sh
    RESULT_VARIABLE RESV
    OUTPUT_VARIABLE OUTV
    ERROR_VARIABLE ERRV)
  #message("res var is ${RESV}")
  #message("out var is ${OUTV}")
  #message("err var is ${ERRV}")

endforeach()

