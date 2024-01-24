#include "catch2/catch.hpp"

#include "physics/spa/spa_functions.hpp"
#include "share/io/scream_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

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
  using SPAFunc = spa::SPAFunctions<Real, DefaultDevice>;

  constexpr Real tol  = std::numeric_limits<Real>::epsilon()*1000;

  // Set up the mpi communicator and init the pio subsystem
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::eam_init_pio_subsystem(comm);

  std::string spa_data_file  = SCREAM_DATA_DIR "/init/spa_data_for_testing.nc";
  std::string spa_remap_file = SCREAM_DATA_DIR "/init/spa_data_for_testing.nc";

  const int ncols_model = 48;
  const int nlevs       = scorpio::get_dimlen(spa_data_file,"lev");
  const int nswbands    = scorpio::get_dimlen(spa_data_file,"swband");
  const int nlwbands    = scorpio::get_dimlen(spa_data_file,"lwband");

  auto grid_model = create_point_grid("model_grid",ncols_model,nlevs,comm);

  // Create horiz remapper
  auto remapper = SPAFunc::create_horiz_remapper(grid_model,spa_data_file,spa_remap_file);
  const int ncols_data = remapper->get_src_grid()->get_num_global_dofs();

  // Create spa data reader
  auto reader = SPAFunc::create_spa_data_reader(remapper,spa_data_file);

  // Recall, SPA data is padded, so we initialize with 2 more levels than the source data file.
  SPAFunc::SPAInput spa_data(grid_model->get_num_local_dofs(), nlevs+2, nswbands, nlwbands);

  // Verify that the interpolated values match the algorithm for the data and the weights.
  //       weights(i) = 1 / (2**i), weights(-1) = 1 / (2**(ncols-1)) such that sum(weights) = 1., for i=0,1,2
  //       FOR t=1,2,3; i=0,1,2; b=1,2 or 1,2,3 and k=0,1,2,3
  //       p(t,i) = (t+1) * (i+1)*100
  //       ccn3(t,i,k) = (i+1)*100 + t*10 + k
  //       aer_g_sw(t,i,b,k) = t
  //       aer_ssa_sw(t,i,b,k) = i
  //       aer_tau_sw(t,i,b,k) = b
  //       aer_tau_lw(t,i,b,k) = k
  auto ps_d         = spa_data.PS;
  auto ccn3_d       = ekat::scalarize(spa_data.data.CCN3);
  auto aer_g_sw_d   = ekat::scalarize(spa_data.data.AER_G_SW);
  auto aer_ssa_sw_d = ekat::scalarize(spa_data.data.AER_SSA_SW);
  auto aer_tau_sw_d = ekat::scalarize(spa_data.data.AER_TAU_SW);
  auto aer_tau_lw_d = ekat::scalarize(spa_data.data.AER_TAU_LW);

  auto ps_h         = Kokkos::create_mirror_view(ps_d);
  auto ccn3_h       = Kokkos::create_mirror_view(ccn3_d);
  auto aer_g_sw_h   = Kokkos::create_mirror_view(aer_g_sw_d);
  auto aer_ssa_sw_h = Kokkos::create_mirror_view(aer_ssa_sw_d);
  auto aer_tau_sw_h = Kokkos::create_mirror_view(aer_tau_sw_d);
  auto aer_tau_lw_h = Kokkos::create_mirror_view(aer_tau_lw_d);

  const int max_time = 3;
  for (int time_index = 0;time_index<max_time; time_index++) {
    SPAFunc::update_spa_data_from_file(*reader, time_index, *remapper, spa_data);

    Kokkos::deep_copy(ps_h,        ps_d);
    Kokkos::deep_copy(ccn3_h,      ccn3_d);
    Kokkos::deep_copy(aer_g_sw_h,  aer_g_sw_d);
    Kokkos::deep_copy(aer_ssa_sw_h,aer_ssa_sw_d);
    Kokkos::deep_copy(aer_tau_sw_h,aer_tau_sw_d);
    Kokkos::deep_copy(aer_tau_lw_h,aer_tau_lw_d);

    for (int idof=0; idof<grid_model->get_num_local_dofs(); ++idof) {
      REQUIRE(std::abs(ps_h(idof) - ps_func(time_index,ncols_data))<tol);
      for (int kk=0; kk<nlevs; kk++) {
        // Recall, SPA data read from file is padded, so we need to offset the kk index for the data by 1.
        REQUIRE(std::abs(ccn3_h(idof,kk+1) - ccn3_func(time_index, kk, ncols_data))<tol);
        for (int n=0; n<nswbands; n++) {
          REQUIRE(aer_g_sw_h(idof,n,kk+1)   == aer_func(time_index,n,kk,ncols_data,0));
          REQUIRE(aer_ssa_sw_h(idof,n,kk+1) == aer_func(time_index,n,kk,ncols_data,1));
          REQUIRE(aer_tau_sw_h(idof,n,kk+1) == aer_func(time_index,n,kk,ncols_data,2));
        }
        for (int n=0; n<nlwbands; n++) {
          REQUIRE(aer_tau_lw_h(idof,n,kk+1) ==  aer_func(time_index,n,kk,ncols_data,3));
        }
      }
    }
  }

  // Clean up
  reader = nullptr;
  scorpio::eam_pio_finalize();
}

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

