set(BETR_MATH_TESTS
  test_betr_mathfunc.pfunit)

configure_file(testSuites.inc testSuites.inc COPYONLY)

set(BETR_MATH_INCLUDE_DIRS
  src/betr/betr_util
  src/betr/betr_math
  )

# turn off compiler warning because we don't have control of pfunit preprocessor
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-use-without-only -Wno-unused-dummy-argument -Wno-error=use-without-only")

include(add_betr_pfunit_test)
add_betr_pfunit_test(test_betr_mathfunc ${BETR_MATH_TESTS} ${BETR_MATH_INCLUDE_DIRS})
