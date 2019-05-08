#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

#include "tridiag_tests.hpp"
#include "share/scream_session.hpp"

int main (int argc, char **argv) {
  int num_failed = 0;
  scream::initialize_scream_session(argc, argv); {
    if (argc > 1) {
      // Performance test.
      scream::tridiag::test::perf::Input in;
      const auto stat = in.parse(argc, argv);
      if (stat) scream::tridiag::test::perf::run<scream::Real>(in);
    } else {
      // Correctness tests.
      num_failed = Catch::Session().run(argc, argv);
    }
  } scream::finalize_scream_session();
  return num_failed != 0 ? 1 : 0;
}
