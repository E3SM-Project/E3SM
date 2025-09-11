#include <catch2/catch.hpp>

#include "share/core/eamxx_pysession.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>

#include <filesystem>

namespace py = pybind11;

namespace scream {

TEST_CASE("pysession", "") {
  auto& ps = PySession::get();

  REQUIRE_THROWS(ps.finalize());

  REQUIRE (not ps.is_initialized());

  ps.initialize();
  REQUIRE (ps.is_initialized());

  REQUIRE_THROWS (ps.safe_import("some_module")); // Can't find module
  ps.add_curr_path();

  // Use a scope, so the py object is destroyed before pysession is finalized
  SECTION ("call_py_fcn") {
    auto my_module = ps.safe_import("my_module");

    auto a = 10;
    auto b = my_module.attr("add_one")(a);
    REQUIRE (b.cast<int>()==(a+1));
  }

  ps.finalize();
  REQUIRE (not ps.is_initialized());
}

} // namespace scream
