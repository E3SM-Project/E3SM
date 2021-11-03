#include "catch2/catch.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "physics/spa/spa_functions.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace {

using namespace scream;
using namespace spa;

template <typename S>
using view_1d = typename KokkosTypes<DefaultDevice>::template view_1d<S>;

TEST_CASE("spa_read_remap_data","spa")
{
  // Set up the mpi communicator and init the pio subsystem
  ekat::Comm spa_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(spa_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Establish the SPA function object
  using SPAFunc         = spa::SPAFunctions<Real, DefaultDevice>;
  SPAFunc::SPAHorizInterp spa_horiz_interp;

  const std::string remap_file_name = "spa_data_for_testing.nc";

  Int tgt_grid_ncols_total = 48;
  Int src_grid_ncols = 20;
  // Break the test set of columns into local degrees of freedom per mpi rank
  auto comm_size = spa_comm.size();
  auto comm_rank = spa_comm.rank();
  std::vector<int> my_dofs;
  for (int ii=comm_rank;ii<tgt_grid_ncols_total;ii+=comm_size) {
    my_dofs.push_back(ii);
  }
  Int tgt_grid_ncols = my_dofs.size();
  // Make sure that the total set of columns has been completely broken up.
  Int test_total_ncols = 0;
  spa_comm.all_reduce(&tgt_grid_ncols,&test_total_ncols,1,MPI_SUM);
  REQUIRE(test_total_ncols == tgt_grid_ncols_total);

  view_1d<int> dofs_gids("",tgt_grid_ncols);
  auto dofs_gids_h = Kokkos::create_mirror_view(dofs_gids);
  Kokkos::deep_copy(dofs_gids_h,dofs_gids);
  for (int ii=0;ii<tgt_grid_ncols;ii++) {
    dofs_gids_h(ii) = my_dofs[ii];
  }
  Kokkos::deep_copy(dofs_gids,dofs_gids_h);
  SPAFunc::get_remap_weights_from_file(remap_file_name,tgt_grid_ncols_total,dofs_gids,spa_horiz_interp);

  REQUIRE(spa_horiz_interp.length==tgt_grid_ncols*src_grid_ncols);
  REQUIRE(spa_horiz_interp.source_grid_ncols==src_grid_ncols);

  // We have a few metrics to ensure that the data read from file matches the data in the file.
  Real tol = 1e5*std::numeric_limits<Real>::epsilon();
  Int col_sum = 0;
  Int row_sum = 0;
  Real wgt_sum = 0.0;
  view_1d<Real> wgts("",tgt_grid_ncols);
  Kokkos::deep_copy(wgts,0.0);
  for (int i=0; i<spa_horiz_interp.length; i++) {
    col_sum += spa_horiz_interp.target_grid_loc[i];
    row_sum += spa_horiz_interp.source_grid_loc[i];
    wgt_sum += spa_horiz_interp.weights[i];
    wgts(spa_horiz_interp.target_grid_loc[i]) += spa_horiz_interp.weights[i];
  }
  // Note, for our test problem the column sum is sum(0...tgt_grid_ncols-1)*src_grid_ncols,
  //       and the                     row sum is sum(0...src_grid_ncols-1)*tgt_grid_ncols. 
  //       Each set of weights should add up to 1.0 for each ncol, so total weights=tgt_grid_ncols
  REQUIRE(col_sum == (tgt_grid_ncols*(tgt_grid_ncols-1))/2*src_grid_ncols);
  REQUIRE(row_sum == (src_grid_ncols*(src_grid_ncols-1)/2*tgt_grid_ncols));
  REQUIRE(std::abs(wgt_sum - 1.0*tgt_grid_ncols) < tol);
  // The sum of remap weights should always be 1.0
  for (int i=0; i<tgt_grid_ncols; i++) {
    REQUIRE(wgts[i]==1.0);
  }
  
  // All Done 
  scorpio::eam_pio_finalize();
}

} //namespace
