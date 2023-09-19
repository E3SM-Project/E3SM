#include <catch2/catch.hpp>
#include "share/io/scream_output_manager.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace {

using namespace scream;

TEST_CASE("create_vert_remap","create_vert_remap")
{
  // Simple function to create a 1D remap column to test nudging w/ remapped data
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  int nlevs = 5*SCREAM_PACK_SIZE+1;
  std::vector<std::int64_t> dofs_levs(nlevs);
  std::iota(dofs_levs.begin(),dofs_levs.end(),0);
  std::vector<Real> p_tgt;
  Real p_top=0, p_bot=102500;
  Real dp = (p_bot - p_top) / (nlevs-1);
  for (int ii=0; ii<nlevs; ++ii) {
    Real p_loc = p_top + dp*ii;
    p_tgt.push_back(p_loc);
  }  

  std::string remap_filename = "vertical_remap.nc";

  scorpio::register_file(remap_filename, scorpio::FileMode::Write);
  scorpio::register_dimension(remap_filename,"lev", "lev", nlevs, false);
  scorpio::register_variable(remap_filename,"p_levs","p_levs","none",{"lev"},"real","real","Real-lev");
  scorpio::set_dof(remap_filename,"p_levs",dofs_levs.size(),dofs_levs.data()); 
  scorpio::eam_pio_enddef(remap_filename);
  scorpio::grid_write_data_array(remap_filename,"p_levs",p_tgt.data(),nlevs);
  scorpio::eam_pio_closefile(remap_filename);
}
} // end namespace
