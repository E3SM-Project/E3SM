// This is a tiny program that calls p3_init() to generate tables used by p3

#include "physics/p3/p3_functions.hpp"
#include "share/core/eamxx_session.hpp"

int main(int argc, char** argv) {
  using P3F = scream::p3::Functions<scream::Real, ekat::DefaultDevice>;

  MPI_Init(&argc, &argv);
  scream::initialize_eamxx_session(argc, argv);
  {
    ekat::Comm comm(MPI_COMM_WORLD);
    P3F::p3_init(/* write_tables = */ true, &comm);
  }
  scream::finalize_eamxx_session();
  MPI_Finalize();

  return 0;
}
