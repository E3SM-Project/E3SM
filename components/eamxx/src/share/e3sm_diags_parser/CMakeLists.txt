cmake_minimum_required(VERSION 3.20)

project(e3sm_diags_parser LANGUAGES CXX)

add_library(e3sm_diags_parser 
src/lexer.cpp
src/ast_print.cpp
src/tokens.cpp
src/precedences.cpp
src/parser.cpp
)

add_executable(edp 
tools/edp.cpp
)
target_link_libraries(edp
    PRIVATE
    e3sm_diags_parser
)

target_compile_features(e3sm_diags_parser PUBLIC cxx_std_20)
target_include_directories(
    e3sm_diags_parser
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
include(CTest)
if(BUILD_TESTING)
    include(FetchContent)

    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.8.1
    )

    FetchContent_MakeAvailable(Catch2)

    add_executable(edp_tests
        tests/test_lexer.cpp
        tests/test_parser.cpp
    )

    target_link_libraries(
        edp_tests
        PRIVATE
            e3sm_diags_parser
            Catch2::Catch2WithMain
    )

    include(Catch)

    catch_discover_tests(edp_tests)
endif()

