set(BETR_CORE_TESTS
  test_betr_transport.pfunit)

# NOTE(bja, 201603) you must add tests to testSuites.inc as
# '${MODULE_NAME}_suite' inorder to register the tests and have them
# run.

configure_file(testSuites.inc testSuites.inc COPYONLY)

set(BETR_CORE_INCLUDE_DIRS
  src/betr/betr_util
  src/betr/betr_core
  )

# turn off compiler warning because we don't have control of pfunit preprocessor
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-use-without-only -Wno-unused-dummy-argument")

include(add_betr_pfunit_test)
add_betr_pfunit_test(test_betr_core ${BETR_CORE_TESTS} ${BETR_CORE_INCLUDE_DIRS})
