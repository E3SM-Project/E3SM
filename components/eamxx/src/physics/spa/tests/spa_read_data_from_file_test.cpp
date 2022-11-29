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

Real ps_func(const int t, const int ncols);
Real ccn3_func(const int t, const int klev, const int ncols);
Real aer_func(const int t, const int bnd, const int klev, const int ncols, const int mode);

TEST_CASE("spa_read_data","spa")
{
  // Set up the mpi communicator and init the pio subsystem
  ekat::Comm spa_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(spa_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Establish the SPA function object
  using SPAFunc = spa::SPAFunctions<Real, DefaultDevice>;
  using Spack = SPAFunc::Spack;
  using gid_type = SPAFunc::gid_type;
  constexpr Real tol  = std::numeric_limits<Real>::epsilon()*1000;

  std::string spa_data_file  = SCREAM_DATA_DIR "/init/spa_data_for_testing.nc";
  std::string spa_remap_file = SCREAM_DATA_DIR "/init/spa_data_for_testing.nc";
  int max_time = 3;
  int ncols    = 48;
  int ncols_src = 20;
  int nlevs    = 4;
  int nswbands = 2;
  int nlwbands = 3;
  // Break the test set of columns into local degrees of freedom per mpi rank
  auto comm_size = spa_comm.size();
  auto comm_rank = spa_comm.rank();

  int my_ncols = ncols/comm_size + (comm_rank < ncols%comm_size ? 1 : 0);
  view_1d<gid_type> dofs_gids("",my_ncols);
  gid_type min_dof = 1; // Start global-ids from 1
  Kokkos::parallel_for("", my_ncols, KOKKOS_LAMBDA(const int& ii) {
    dofs_gids(ii) = min_dof + static_cast<gid_type>(comm_rank + ii*comm_size);
  });
  // Make sure that the total set of columns has been completely broken up.
  int test_total_ncols = 0;
  spa_comm.all_reduce(&my_ncols,&test_total_ncols,1,MPI_SUM);
  REQUIRE(test_total_ncols == ncols);

  // Set up the set of SPA structures needed to run the test
  SPAFunc::SPAHorizInterp spa_horiz_interp;
  spa_horiz_interp.m_comm = spa_comm;
  SPAFunc::get_remap_weights_from_file(spa_remap_file,min_dof,dofs_gids,spa_horiz_interp);
  // Recall, SPA data is padded, so we initialize with 2 more levels than the source data file.
  SPAFunc::SPAInput spa_data(dofs_gids.size(), nlevs+2, nswbands, nlwbands);

  // Verify that the interpolated values match the algorithm for the data and the weights.
  //       weights(i) = 1 / (2**i), weights(-1) = 1 / (2**(ncols-1)) such that sum(weights) = 1., for i=0,1,2
  //       FOR t=1,2,3; i=0,1,2; b=1,2 or 1,2,3 and k=0,1,2,3
  //       p(t,i) = (t+1) * (i+1)*100
  //       ccn3(t,i,k) = (i+1)*100 + t*10 + k
  //       aer_g_sw(t,i,b,k) = t
  //       aer_ssa_sw(t,i,b,k) = i
  //       aer_tau_sw(t,i,b,k) = b
  //       aer_tau_lw(t,i,b,k) = k
  auto ps_h         = Kokkos::create_mirror_view(spa_data.PS);
  auto ccn3_h       = Kokkos::create_mirror_view(spa_data.data.CCN3);
  auto aer_g_sw_h   = Kokkos::create_mirror_view(spa_data.data.AER_G_SW);
  auto aer_ssa_sw_h = Kokkos::create_mirror_view(spa_data.data.AER_SSA_SW);
  auto aer_tau_sw_h = Kokkos::create_mirror_view(spa_data.data.AER_TAU_SW);
  auto aer_tau_lw_h = Kokkos::create_mirror_view(spa_data.data.AER_TAU_LW);
  auto dofs_gids_h = Kokkos::create_mirror_view(dofs_gids);
  Kokkos::deep_copy(dofs_gids_h,dofs_gids);
  for (int time_index = 0;time_index<max_time; time_index++) {
    SPAFunc::update_spa_data_from_file(spa_data_file, time_index, nswbands, nlwbands,
                                       spa_horiz_interp, spa_data);
    Kokkos::deep_copy(ps_h,spa_data.PS);
    Kokkos::deep_copy(ccn3_h,spa_data.data.CCN3);
    Kokkos::deep_copy(aer_g_sw_h,spa_data.data.AER_G_SW);
    Kokkos::deep_copy(aer_ssa_sw_h,spa_data.data.AER_SSA_SW);
    Kokkos::deep_copy(aer_tau_sw_h,spa_data.data.AER_TAU_SW);
    Kokkos::deep_copy(aer_tau_lw_h,spa_data.data.AER_TAU_LW);
    for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
      REQUIRE(std::abs(ps_h(dof_i) - ps_func(time_index,ncols_src))<tol);
      for (int kk=0;kk<nlevs;kk++) {
        // Recall, SPA data read from file is padded, so we need to offset the kk index for the data by 1.
        int kpack = (kk+1) / Spack::n;
        int kidx  = (kk+1) % Spack::n;
        REQUIRE(std::abs(ccn3_h(dof_i,kpack)[kidx] - ccn3_func(time_index, kk, ncols_src))<tol);
        for (int n=0;n<nswbands;n++) {
          REQUIRE(aer_g_sw_h(dof_i,n,kpack)[kidx]   == aer_func(time_index,n,kk,ncols_src,0));
          REQUIRE(aer_ssa_sw_h(dof_i,n,kpack)[kidx] == aer_func(time_index,n,kk,ncols_src,1));
          REQUIRE(aer_tau_sw_h(dof_i,n,kpack)[kidx] == aer_func(time_index,n,kk,ncols_src,2));
        }
        for (int n=0;n<nlwbands;n++) {
          REQUIRE(aer_tau_lw_h(dof_i,n,kpack)[kidx] ==  aer_func(time_index,n,kk,ncols_src,3));
        }
      }
    }
  }
  
  // All Done 
  scorpio::eam_pio_finalize();
} // run_property

// Some helper functions for the require statements:
Real ps_func(const int t, const int ncols)
{
  Real ps = 0.0;
  for (int i=1;i<=ncols;i++) {
    Real wgt = 1.0 / std::pow(2.0,i);
    if (i == ncols) { wgt *= 2.0; }
    ps += (t+1) * i*100.0 * wgt;
  }
  return ps;
} // ps_func
//
Real ccn3_func(const int t, const int klev, const int ncols)
{
  Real ccn3 = 0.0;
  for (int i=1;i<=ncols;i++) {
    Real wgt = 1.0 / std::pow(2.0,i);
    if (i == ncols) {wgt *= 2.0;}
    ccn3 += wgt * (klev*1.0 + t*10.0 + i*100.0);
  }
  return ccn3;
} // ccn3_func
//
Real aer_func(const int t, const int bnd, const int klev, const int ncols, const int mode)
{
  Real aer_out = 0.0;
  for (int i=1;i<=ncols;i++) {
    Real wgt = 1.0 / std::pow(2.0,i);
    if (i == ncols) {wgt *= 2.0;}
    if (mode==0) {  // G
      aer_out += wgt * t;
    } else if (mode==1) { // SSA
      aer_out += wgt * (i-1);
    } else if (mode==2) { // TAU
      aer_out += wgt * bnd;
    } else if (mode==3) { // LW-TAU
      aer_out += wgt * klev;
    }
  }
  return aer_out;
} // aer_func

} // namespace

