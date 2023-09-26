#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include <Kokkos_Core.hpp>

#include <Hommexx_Session.hpp>

#include <Context.hpp>

// Make command-line arguments available to tests. Anything following "hommexx"
// in the command-line argument list is passed to the tests.
int hommexx_catch2_argc;
char** hommexx_catch2_argv;

int main(int argc, char **argv) {

  Homme::initialize_hommexx_session();
  // Filter arguments so catch2 doesn't try to interpret hommexx-specific ones.
  hommexx_catch2_argc = argc;
  hommexx_catch2_argv = argv;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "hommexx") {
      argc = i;
      hommexx_catch2_argc -= i + 1;
      hommexx_catch2_argv = argv + i + 1;
      break;
    }
  }

  int result = Catch::Session().run(argc, argv);

  return result;
}
