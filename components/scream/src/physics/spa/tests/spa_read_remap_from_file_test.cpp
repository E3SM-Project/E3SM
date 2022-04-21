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
template <typename S>
using view_1d_host = typename view_1d<S>::HostMirror; //KokkosTypes<HostDevice>::template view_1d<S>;

TEST_CASE("spa_read_remap_data","spa")
{
  // Set up the mpi communicator and init the pio subsystem
  ekat::Comm spa_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(spa_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Establish the SPA function object
  using SPAFunc  = spa::SPAFunctions<Real, DefaultDevice>;
  using gid_type = SPAFunc::gid_type;
  SPAFunc::SPAHorizInterp spa_horiz_interp;
  spa_horiz_interp.m_comm = spa_comm;

  const std::string remap_file_name = SCREAM_DATA_DIR "/init/spa_data_for_testing.nc";

  Int tgt_grid_ncols_total = 48;
  Int src_grid_ncols = 20;
  // Break the test set of columns into local degrees of freedom per mpi rank
  auto comm_size = spa_comm.size();
  auto comm_rank = spa_comm.rank();
  int tgt_grid_ncols = tgt_grid_ncols_total/comm_size + (comm_rank < tgt_grid_ncols_total%comm_size ? 1 : 0);
  view_1d<gid_type> dofs_gids("",tgt_grid_ncols);
  gid_type min_dof = 0;  // We will set up the dof's to start from 0
  Kokkos::parallel_for("", tgt_grid_ncols, KOKKOS_LAMBDA(const int& ii) {
    dofs_gids(ii) = min_dof + static_cast<gid_type>(comm_rank + ii*comm_size);
  });
  // Make sure that the total set of columns has been completely broken up.
  Int test_total_ncols = 0;
  spa_comm.all_reduce(&tgt_grid_ncols,&test_total_ncols,1,MPI_SUM);
  REQUIRE(test_total_ncols == tgt_grid_ncols_total);

  SPAFunc::get_remap_weights_from_file(remap_file_name,tgt_grid_ncols_total,0,dofs_gids,spa_horiz_interp);
  // Make sure one_to_one remap has the correct unique columns
  REQUIRE(spa_horiz_interp.num_unique_cols==src_grid_ncols);

  REQUIRE(spa_horiz_interp.length==tgt_grid_ncols*src_grid_ncols);
  REQUIRE(spa_horiz_interp.source_grid_ncols==src_grid_ncols);

  // We have a few metrics to ensure that the data read from file matches the data in the file.
  Int col_sum = 0;
  Int row_sum = 0;
  view_1d_host<Real> wgts("",tgt_grid_ncols);
  Kokkos::deep_copy(wgts,0.0);
  auto weights_h = Kokkos::create_mirror_view(spa_horiz_interp.weights);
  Kokkos::deep_copy(weights_h,spa_horiz_interp.weights);
  auto target_h = Kokkos::create_mirror_view(spa_horiz_interp.target_grid_loc);
  Kokkos::deep_copy(target_h,spa_horiz_interp.target_grid_loc);
  for (int i=0; i<spa_horiz_interp.length; i++) {
    wgts(target_h[i]) += weights_h[i];
  }
  Kokkos::parallel_reduce("", spa_horiz_interp.length, KOKKOS_LAMBDA (const int& i, Int& lsum) {
    lsum += spa_horiz_interp.target_grid_loc[i];
  },row_sum);
  Kokkos::parallel_reduce("", spa_horiz_interp.length, KOKKOS_LAMBDA (const int& i, Int& lsum) {
    lsum += spa_horiz_interp.source_grid_loc[i];
  },col_sum);
  // Note, for our test problem the column sum is sum(0...tgt_grid_ncols-1)*src_grid_ncols,
  //       and the                     row sum is sum(0...src_grid_ncols-1)*tgt_grid_ncols. 
  //       Each set of weights should add up to 1.0 for each ncol, so total weights=tgt_grid_ncols
  REQUIRE(row_sum == (tgt_grid_ncols*(tgt_grid_ncols-1))/2*src_grid_ncols);
  REQUIRE(col_sum == (src_grid_ncols*(src_grid_ncols-1)/2*tgt_grid_ncols));
  // The sum of remap weights should always be 1.0
  Real wgt_sum = 0.0;
  for (int i=0; i<tgt_grid_ncols; i++) {
    REQUIRE(wgts[i]==1.0);
    wgt_sum += wgts[i];
  }
  REQUIRE(wgt_sum/tgt_grid_ncols == 1.0);
  
  // All Done 
  scorpio::eam_pio_finalize();
}

} //namespace
