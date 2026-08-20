# Prefer a Catch2 that is already installed; fall back to fetching one, so that
# a bare `cmake -S . -B build` still works on a machine that has never seen it.
find_package(Catch2 3 QUIET)

if (NOT Catch2_FOUND)
  include(FetchContent)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
  )

  FetchContent_MakeAvailable(Catch2)

  # Needed for include(Catch) below when Catch2 comes in as a subproject.
  list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
endif()

add_executable(edp_tests
  test_lexer.cpp
  test_parser.cpp
)

target_link_libraries(edp_tests
  PRIVATE
    e3sm_diags_parser
    Catch2::Catch2WithMain
)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(edp_tests PRIVATE -Wall -Wextra)
endif()

include(Catch)

# Registers each TEST_CASE as its own ctest test, so a failure names itself.
catch_discover_tests(edp_tests)
