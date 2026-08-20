cmake_minimum_required(VERSION 3.20)

project(e3sm_diags_parser
  VERSION 0.1.0
  DESCRIPTION "Lexer and parser for the diagnostics expression language"
  LANGUAGES CXX)

# This library is deliberately standalone: it has no Kokkos, MPI, netCDF or
# EKAT dependency, and is not wired into the CIME/csm_share build. That is
# what lets it be configured and tested anywhere cmake and a C++20 compiler
# exist, including a stock GitHub runner.
option(EDP_ENABLE_TESTS "Build the unit tests" ON)
option(EDP_ENABLE_TOOL  "Build the command line tool" ON)

add_library(e3sm_diags_parser
  src/ast_print.cpp
  src/lexer.cpp
  src/parser.cpp
  src/precedences.cpp
  src/tokens.cpp
)

target_compile_features(e3sm_diags_parser PUBLIC cxx_std_20)

target_include_directories(e3sm_diags_parser
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Warnings stay on: with no external headers in the include path there is
# nothing here that anyone has to fight, so the library should stay clean.
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(e3sm_diags_parser PRIVATE -Wall -Wextra -Wpedantic)
endif()

if (EDP_ENABLE_TOOL)
  add_executable(edp tools/edp.cpp)
  target_link_libraries(edp PRIVATE e3sm_diags_parser)
endif()

if (EDP_ENABLE_TESTS)
  include(CTest)
  if (BUILD_TESTING)
    add_subdirectory(tests)
  endif()
endif()
