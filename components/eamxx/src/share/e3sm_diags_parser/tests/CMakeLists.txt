set (EDP_TEST_SRCS
  test_lexer.cpp
  test_parser.cpp
)

if (EDP_STANDALONE)
  # Nothing else provides Catch2 in a standalone build, so fetch one.
  include(FetchContent)

  FetchContent_Declare(
      Catch2
      GIT_REPOSITORY https://github.com/catchorg/Catch2.git
      GIT_TAG v3.8.1
  )

  FetchContent_MakeAvailable(Catch2)

  add_executable(edp_tests ${EDP_TEST_SRCS})

  target_link_libraries(
      edp_tests
      PRIVATE
          e3sm_diags_parser
          Catch2::Catch2WithMain
  )

  include(Catch)

  catch_discover_tests(edp_tests)

# NOTE: tests inside this if statement won't be built in a baselines-only build
elseif (NOT SCREAM_ONLY_GENERATE_BASELINES)
  # Deliberately NOT using CreateUnitTest here. That would link eamxx_test_support,
  # dragging in eamxx_core -> ekat -> Kokkos, so building this folder alone would
  # build Kokkos first. The parser is plain C++ with no Kokkos and no MPI, so we
  # link the (header-only) Catch2 that EKAT provides, supply our own main, and
  # register the test directly. It still runs under ctest, and hence under
  # test-all-eamxx, like every other unit test.
  add_executable(edp_tests ${EDP_TEST_SRCS} edp_catch_main.cpp)

  target_link_libraries(
      edp_tests
      PRIVATE
          e3sm_diags_parser
          Catch2::Catch2
  )

  add_test(NAME edp_tests COMMAND edp_tests)
  set_tests_properties(edp_tests PROPERTIES LABELS parser PROCESSORS 1)
endif()
