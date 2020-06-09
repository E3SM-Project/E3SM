include(EkatCreateUnitTest)

# Test packs
set (PACK_TESTS_SRCS
  pack_tests.cpp
  pack_kokkos_tests.cpp
)
EkatCreateUnitTest(packs "${PACK_TESTS_SRCS}"
  LIBS ekat
  INCLUDE_DIRS "${EKAT_INCLUDE_DIRS}"
)
