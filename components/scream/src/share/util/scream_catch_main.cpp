#define CATCH_CONFIG_RUNNER

#include "catch2/catch.hpp"

#include "share/scream_session.hpp"

int main (int argc, char **argv) {
  int result = -1;
  scream::initialize_scream_session(argc, argv); {
    result = Catch::Session().run(argc, argv);
  } scream::finalize_scream_session();
  return result;
}
