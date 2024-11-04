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

using SPAFunc = spa::SPAFunctions<Real, DefaultDevice>;
using gid_type = SPAFunc::gid_type;

Real ps_func  (const int t, const gid_type icol);
Real ccn3_func(const int t, const gid_type icol, const int klev);
Real aer_func (const int t, const gid_type icol, const int bnd, const int klev, const int mode);

TEST_CASE("spa_one_to_one_remap","spa")
{
  // Set up the mpi communicator and init the pio subsystem
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  std::string spa_data_file = SCREAM_DATA_DIR "/init/spa_data_for_testing.nc";

  const int ncols_model = scorpio::get_dimlen(spa_data_file,"ncol");
  const int nlevs       = scorpio::get_dimlen(spa_data_file,"lev");
  const int nswbands    = scorpio::get_dimlen(spa_data_file,"swband");
  const int nlwbands    = scorpio::get_dimlen(spa_data_file,"lwband");

  auto grid_model = create_point_grid("model_grid",ncols_model,nlevs,comm);

  // Create horiz remapper
  auto remapper = SPAFunc::create_horiz_remapper(grid_model,spa_data_file,"SHOULD_NOT_BE_NEEDED");
  const int ncols_data = remapper->get_src_grid()->get_num_global_dofs();

  // This test should be tailored to have same number of columns as the spa data
  REQUIRE (ncols_data==ncols_model);
  REQUIRE (std::dynamic_pointer_cast<IdentityRemapper>(remapper));

  // Create spa data reader
  auto reader = SPAFunc::create_spa_data_reader(remapper,spa_data_file);
  // TODO: We can add lat/lon to spa_data_file and test the impl for IOP
  util::TimeStamp dummy_ts;
  std::shared_ptr<SPAFunc::IOPReader> dummy_iop_reader;

  // Recall, SPA data is padded, so we initialize with 2 more levels than the source data file.
  SPAFunc::SPAInput spa_data(grid_model->get_num_local_dofs(), nlevs+2, nswbands, nlwbands);

  // Verify that the interpolated values match the algorithm for the data and the weights.
  //       weights(i) = 1.0
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

  auto dofs_gids_h  = grid_model->get_dofs_gids().get_view<const gid_type*,Host>();

  const int max_time = 3;
  for (int time_index = 0;time_index<max_time; time_index++) {
    SPAFunc::update_spa_data_from_file(reader, dummy_iop_reader, dummy_ts, time_index, *remapper, spa_data);
    Kokkos::deep_copy(ps_h,ps_d);
    Kokkos::deep_copy(ccn3_h,ccn3_d);
    Kokkos::deep_copy(aer_g_sw_h,aer_g_sw_d);
    Kokkos::deep_copy(aer_ssa_sw_h,aer_ssa_sw_d);
    Kokkos::deep_copy(aer_tau_sw_h,aer_tau_sw_d);
    Kokkos::deep_copy(aer_tau_lw_h,aer_tau_lw_d);
    for (int idof=0; idof<grid_model->get_num_local_dofs(); ++idof) {
      gid_type glob_i = dofs_gids_h(idof);
      REQUIRE(ps_h(idof) == ps_func(time_index,glob_i));
      for (int kk=0; kk<nlevs; kk++) {
        // Recall, SPA data read from file is padded, so we need to offset the kk index for the data by 1.
        REQUIRE(ccn3_h(idof,kk+1) == ccn3_func(time_index, glob_i, kk));
        for (int n=0; n<nswbands; n++) {
          REQUIRE(aer_g_sw_h(idof,n,kk+1)   == aer_func(time_index,glob_i,n,kk,0));
          REQUIRE(aer_ssa_sw_h(idof,n,kk+1) == aer_func(time_index,glob_i,n,kk,1));
          REQUIRE(aer_tau_sw_h(idof,n,kk+1) == aer_func(time_index,glob_i,n,kk,2));
        }
        for (int n=0; n<nlwbands; n++) {
          REQUIRE(aer_tau_lw_h(idof,n,kk+1) ==  aer_func(time_index,glob_i,n,kk,3));
        }
      }
    }
  }

  // All Done
  reader = nullptr;
  scorpio::finalize_subsystem();
} // run_property

// Some helper functions for the require statements:
Real ps_func(const int t, const gid_type i)
{
  return (t+1) * (i+1)*100.0;
} // ps_func
//
Real ccn3_func(const int t, const gid_type i, const int klev)
{
  return (klev*1.0 + t*10.0 + (i+1)*100.0);
} // ccn3_func
//
Real aer_func(const int t, const gid_type i, const int bnd, const int klev, const int mode)
{
  if (mode==0) {  // G
    return  t;
  } else if (mode==1) { // SSA
    return  i;
  } else if (mode==2) { // TAU
    return  bnd;
  } else if (mode==3) { // LW-TAU
    return  klev;
  }
  EKAT_REQUIRE_MSG(false,"Error in aer_func, incorrect mode passed"); // Shouldn't get here, but just in case
} // aer_func

} // namespace
