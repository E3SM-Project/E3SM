// This is a tiny program that calls p3_init() to generate tables used by p3

#include "physics/p3/p3_functions.hpp"
#include "share/scream_session.hpp"

int main(int argc, char** argv) {
  using P3F = scream::p3::Functions<scream::Real, ekat::DefaultDevice>;

  scream::initialize_scream_session(argc, argv);
  P3F::p3_init(/* write_tables = */ true, /* masterproc */ true);
  scream::finalize_scream_session();

  return 0;
}