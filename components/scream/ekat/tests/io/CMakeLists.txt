include(EkatCreateUnitTest)

# Test utilities (f90)
set (ARRAY_IO_SRCS
  array_io_tests.cpp
  unit_test.f90
)
EkatCreateUnitTest(array_io "${ARRAY_IO_SRCS}"
  LIBS ekat
  INCLUDE_DIRS ${EKAT_INCLUDE_DIRS}
)

# YAML parser
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.yaml ${CMAKE_CURRENT_BINARY_DIR}/input.yaml COPYONLY)
EkatCreateUnitTest(yaml_parser yaml_parser.cpp
  LIBS ekat
  INCLUDE_DIRS ${EKAT_INCLUDE_DIRS}
)
