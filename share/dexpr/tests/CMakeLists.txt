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

add_executable(dexpr_tests
  test_function_registry.cpp
  test_lexer.cpp
  test_parser.cpp
)

target_link_libraries(dexpr_tests
  PRIVATE
    dexpr
    Catch2::Catch2WithMain
)

target_compile_options(dexpr_tests PRIVATE ${DEXPR_WARNINGS})

include(Catch)

# Registers each TEST_CASE as its own ctest test, so a failure names itself.
catch_discover_tests(dexpr_tests)
