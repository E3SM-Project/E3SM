# The diagnostics DSL parser can be built in two ways:
#
#  - standalone, driven by run_tests.sh, where this is the top-level CMake
#    project and fetches its own copy of Catch2;
#  - as part of EAMxx, where src/share/CMakeLists.txt pulls this folder in and
#    the tests reuse the Catch2/test-session machinery shared by all other
#    EAMxx unit tests (so they show up in ctest, and hence test-all-eamxx).
#
# This check must happen before the project() call below.
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set (EDP_STANDALONE TRUE)
else()
  set (EDP_STANDALONE FALSE)
endif()

if (EDP_STANDALONE)
  cmake_minimum_required(VERSION 3.20)
  project(e3sm_diags_parser LANGUAGES CXX)
endif()

add_library(e3sm_diags_parser
src/lexer.cpp
src/ast_print.cpp
src/tokens.cpp
src/precedences.cpp
src/parser.cpp
)

target_compile_features(e3sm_diags_parser PUBLIC cxx_std_20)
target_include_directories(
    e3sm_diags_parser
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# The 'edp' tool is a developer convenience; skip it in library-only (CIME) builds.
if (EDP_STANDALONE OR NOT SCREAM_LIB_ONLY)
  add_executable(edp
  tools/edp.cpp
  )
  target_link_libraries(edp
      PRIVATE
      e3sm_diags_parser
  )
endif()

if (EDP_STANDALONE)
  include(CTest)
  if (BUILD_TESTING)
    add_subdirectory(tests)
  endif()
elseif (NOT SCREAM_LIB_ONLY)
  add_subdirectory(tests)
endif()
