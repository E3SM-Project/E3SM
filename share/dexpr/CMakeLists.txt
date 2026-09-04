cmake_minimum_required(VERSION 3.14)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(dexpr
    VERSION 0.1.0
    DESCRIPTION "Lexer and parser for the diagnostics expression language"
    LANGUAGES CXX)
endif()

# This library is deliberately standalone: it has no Kokkos, MPI, netCDF or
# EKAT dependency, and is not wired into the CIME/csm_share build. That is
# what lets it be configured and tested anywhere cmake and a C++20 compiler
# exist, including a stock GitHub runner.
option(DEXPR_ENABLE_TESTS "Build the unit tests" ON)
option(DEXPR_ENABLE_TOOL  "Build the command line tool" ON)

add_library(dexpr
  src/ast_print.cpp
  src/lexer.cpp
  src/parser.cpp
  src/precedences.cpp
  src/supported_functions.cpp
  src/tokens.cpp
)

target_compile_features(dexpr PUBLIC cxx_std_20)

target_include_directories(dexpr
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

option(DEXPR_WERROR "Treat compiler warnings as errors" OFF)

# Warnings stay on: with no external headers in the include path there is
# nothing here that anyone has to fight, so the library should stay clean.
# CI turns DEXPR_WERROR on to keep it that way. The flags are set here rather
# than per target so tests/ inherits exactly the same set.
set(DEXPR_WARNINGS)
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  set(DEXPR_WARNINGS -Wall -Wextra -Wpedantic)
  if (DEXPR_WERROR)
    list(APPEND DEXPR_WARNINGS -Werror)
  endif()
endif()

target_compile_options(dexpr PRIVATE ${DEXPR_WARNINGS})

if (DEXPR_ENABLE_TOOL)
  # The target is dexpr_cli because the library already owns the name dexpr;
  # the binary it produces is still just `dexpr`.
  add_executable(dexpr_cli tools/dexpr.cpp)
  set_target_properties(dexpr_cli PROPERTIES OUTPUT_NAME dexpr)
  target_link_libraries(dexpr_cli PRIVATE dexpr)
endif()

if (DEXPR_ENABLE_TESTS)
  include(CTest)
  if (BUILD_TESTING)
    add_subdirectory(tests)
  endif()
endif()
