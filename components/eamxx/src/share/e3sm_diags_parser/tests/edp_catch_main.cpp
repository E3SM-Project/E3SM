// The Catch2 that EKAT provides is header-only, so the test executable has to
// define its own entry point. We deliberately do not reuse EKAT's catch main:
// it initializes MPI and Kokkos, which the parser has no use for.
//
// Standalone builds link Catch2::Catch2WithMain instead and never compile this
// file (see tests/CMakeLists.txt).
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
