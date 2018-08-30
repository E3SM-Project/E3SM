#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include "scream_kokkos.hpp"
#include "scream_util.hpp"

int main(int argc, char **argv) {
  int result = -1;
  scream::util::initialize(argc, argv);
  { result = Catch::Session().run(argc, argv); }
  scream::util::finalize();
  return result;
}
