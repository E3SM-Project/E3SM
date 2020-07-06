#include <catch2/catch.hpp>

#include "ekat/scream_pack.hpp"
#include "mct_coupling/scream_scorpio_interface.hpp"
#include "ekat/mpi/scream_comm.hpp"
#include "ekat/scream_types.hpp"


namespace {

TEST_CASE("scorpio_interface", "") {

  using namespace scream;
  using namespace scream::scorpio;

  // Create the set of SCORPIO output files and their respective
  // dimensions and variables.
  int compid=0;
  MPI_Fint fcomm = MPI_Comm_c2f(MPI_COMM_WORLD);
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  // Register the set of output files:
  std::string outfilename = "scorpio_output_test.nc";
  register_outfile(outfilename);
  // Register the set of dimensions per output file
  int xlen=10, ylen=3, zlen=2;
  register_dimension(outfilename,"x","horizontal distance",xlen);
  register_dimension(outfilename,"y","vertical distance",ylen);
  register_dimension(outfilename,"z","height",zlen);
  register_dimension(outfilename,"time","time",0);
  // Finished with the initialization of variables in output file
  eam_pio_enddef();

}

} //namespace
