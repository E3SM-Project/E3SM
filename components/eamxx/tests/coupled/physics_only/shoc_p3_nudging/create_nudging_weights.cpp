#include <catch2/catch.hpp>
#include "share/io/scream_output_manager.hpp"
#include "share/io/scream_scorpio_interface.hpp"

namespace {

using namespace scream;

void create_nudging_weights_ncfile(int ncols, int nlevs, const std::string& filename)
{
  // Simple function to create a 1D remap column to test nudging w/ remapped data
  ekat::Comm io_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(io_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  std::vector<std::int64_t> dofs(ncols*nlevs);
  std::iota(dofs.begin(),dofs.end(),0);
  Real plev[nlevs];
  Real p_top=0, p_bot=102500;
  Real dp = (p_bot - p_top) / (nlevs-1);
  for (int ii=0; ii<nlevs; ++ii) {
    plev[ii] = p_top + dp*ii;
  }  

  Real weights[nlevs][ncols];
  for (auto ilev=0; ilev<nlevs; ++ilev) {
      for (auto icol=0; icol<ncols; ++icol) {
	  if (plev[ilev] <= 1.0e5 && plev[ilev] >= 8.0e4) {
	     weights[ilev][icol] = 1.;
	  } else {
	     weights[ilev][icol] = 0.;
	  }
       }
   }

  scorpio::register_file(filename, scorpio::FileMode::Write);
  scorpio::register_dimension(filename,"ncol", "ncol", ncols, false);
  scorpio::register_dimension(filename,"nlev", "nlev", nlevs, false);
  scorpio::register_variable(filename,"nudging_weights","nudging_weights","none",{"nlev","ncol"},"real","real","Real-lev");
  scorpio::set_dof(filename,"nudging_weights",dofs.size(),dofs.data());
  scorpio::eam_pio_enddef(filename);
  scorpio::grid_write_data_array(filename,"nudging_weights",weights[0],ncols*nlevs);
  scorpio::eam_pio_closefile(filename);
  scorpio::eam_pio_finalize();
}

TEST_CASE("create_nudging_weights","create_nudging_weights")
{
  create_nudging_weights_ncfile(218, 6, "nudging_weights_remapped.nc");
  create_nudging_weights_ncfile(218, 128, "nudging_weights.nc");
}
} // end namespace
