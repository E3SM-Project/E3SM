#include <catch2/catch.hpp>
#include "share/io/eamxx_output_manager.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

namespace {

using namespace scream;

void create_vert_remap() {
  // Simple function to create a 1D remap column to test nudging w/ remapped data
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  int nlevs = 5*SCREAM_PACK_SIZE+1;
  std::vector<Real> p_tgt;
  Real p_top=0, p_bot=102500;
  Real dp = (p_bot - p_top) / (nlevs-1);
  for (int ii=0; ii<nlevs; ++ii) {
    Real p_loc = p_top + dp*ii;
    p_tgt.push_back(p_loc);
  }

  std::string remap_filename = "vertical_remap.nc";

  scorpio::register_file(remap_filename, scorpio::FileMode::Write);
  scorpio::define_dim(remap_filename, "lev", nlevs);
  scorpio::define_var(remap_filename,"p_levs",{"lev"},"real");
  scorpio::enddef(remap_filename);
  scorpio::write_var(remap_filename,"p_levs",p_tgt.data());
  scorpio::release_file(remap_filename);
  scorpio::finalize_subsystem();
}

void create_nudging_weights_ncfile(int ntimes, int ncols, int nlevs, const std::string& filename)
{
  // Simple function to create a 1D remap column to test nudging w/ remapped data
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  Real plev[nlevs];
  Real p_top=0, p_bot=102500;
  Real dp = (p_bot - p_top) / (nlevs-1);
  for (int ilev=0; ilev<nlevs; ++ilev) {
    plev[ilev] = p_top + dp*ilev;
  }  

  Real weights[ntimes][ncols][nlevs];
  for (auto itime=0; itime<ntimes; ++itime) {
    for (auto ilev=0; ilev<nlevs; ++ilev) {
       for (auto icol=0; icol<ncols; ++icol) {
	  if (plev[ilev] <= 1.0e5 && plev[ilev] >= 8.0e4) {
	     weights[itime][icol][ilev] = 1.;
	  } else {
	     weights[itime][icol][ilev] = 0.;
	  }
       }
    }
  }

  // NOTE: we are generating a file with "time" dimension that is not unlimited.
  //       This should stress test AtmosphereInput, which will 'pretend' that
  //       a dimension called 'time' is unlimited, even though it isn't in the file
  scorpio::register_file(filename, scorpio::FileMode::Write);
  scorpio::define_dim(filename,"ncol",ncols);
  scorpio::define_dim(filename,"lev", nlevs);
  scorpio::define_dim(filename,"time",ntimes);
  scorpio::define_var(filename,"nudging_weights",{"time","ncol","lev"},"real");
  scorpio::enddef(filename);
  scorpio::write_var(filename,"nudging_weights",&weights[0][0][0]);
  scorpio::release_file(filename);
  scorpio::finalize_subsystem();
}

TEST_CASE("create_vert_remap_and_weights","create_vert_remap_and_weights")
{
  create_vert_remap();
  create_nudging_weights_ncfile(1, 218, 72, "nudging_weights.nc");
  create_nudging_weights_ncfile(1, 218, 128, "nudging_weights_L128.nc");
}
} // end namespace
