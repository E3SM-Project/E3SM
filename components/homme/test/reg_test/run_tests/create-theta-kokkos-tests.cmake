##### this is needed if we want dereferencing work as now in test-lists.sh
#CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#theta-form0-ttype5-hvs1-hvst0-r3-q1-nutop0-samenu-template


SET(BUILD_HOMME_THETA_KOKKOS TRUE)
message("${BUILD_HOMME_THETA_KOKKOS}")
include(test-list.cmake)


message("${HOMME_THETA_TESTS_WITH_PROFILE_1}")

foreach(jj ${HOMME_THETA_TESTS_WITH_PROFILE_1})

message(${jj})

#parse jj
string(REPLACE "-" ";" jjlist ${jj})

message("${jjlist}")

#now we can access jj encoded stuff, samenu and nutop a exceptions
#theta;form0;ttype5;hvs1;hvst0;r3;q1;nutop0;samenu

#form:
  list(GET jjlist 1 aa)
#now get a value from aa
  string(REGEX MATCH "[0-9]" bb "${aa}")
  message("adv form is ${bb}")
  set(AAFORM "${bb}")

#ttype:
  list(GET jjlist 2 aa)
#now get a value from aa
  string(REGEX MATCH "[0-9]" bb "${aa}")
  message("ttype is ${bb}")
  set(AATTYPE "${bb}")

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
  set(AATSVH "${bb}")

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

set(AANAME ${jj})
configure_file(mmm.sh mmm2.sh)



#now will all that configure template
#file (READ "theta-xx-test-template.cmake" ffile)
#file(REPLACE)


#  string(REGEX MATCHALL "" llist "${NAME}")
#  configure_file(${jj}-template.cmake ${jj}.cmake COPYONLY)
#list()

endforeach()

