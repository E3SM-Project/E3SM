#ifndef EDP_TESTS_EDP_CATCH_HPP
#define EDP_TESTS_EDP_CATCH_HPP

// These tests run against two different Catch2 releases: v3 when the parser is
// built standalone (fetched by tests/CMakeLists.txt) and v2 when it is built as
// part of EAMxx (provided by EKAT). The macros we use are the same in both, but
// the headers that declare them are not, so pick whichever is available.
#if __has_include(<catch2/catch_test_macros.hpp>)
#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>
#else
#include <catch2/catch.hpp>
#endif

#endif // EDP_TESTS_EDP_CATCH_HPP
