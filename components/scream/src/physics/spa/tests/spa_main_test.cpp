#include "catch2/catch.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "physics/spa/spa_functions.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

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
using view_2d = typename KokkosTypes<DefaultDevice>::template view_2d<S>;
template <typename S>
using view_3d = typename KokkosTypes<DefaultDevice>::template view_3d<S>;

using SPAFunc = spa::SPAFunctions<Real, DefaultDevice>;
using Spack = SPAFunc::Spack;

// Helper Functions
void compute_max_min(const view_1d<const Spack>& input, const int start, const int end, Real& min, Real& max);

TEST_CASE("spa_read_data","spa")
{
  using C = scream::physics::Constants<Real>;
  static constexpr auto P0 = C::P0;

 /* ====================================================================
  *                  Test Setup, create structures, etc.
  * ====================================================================*/

  // Set up the mpi communicator and init the pio subsystem
  ekat::Comm spa_comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  MPI_Fint fcomm = MPI_Comm_c2f(spa_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  scorpio::eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler

  // Establish the SPA function object
  using SPAFunc = spa::SPAFunctions<Real, DefaultDevice>;
  using Spack = SPAFunc::Spack;
  using gid_type = SPAFunc::gid_type;

  std::string fname = "spa_main.yaml";
  ekat::ParameterList test_params("Atmosphere Driver");
  REQUIRE_NOTHROW ( parse_yaml_file(fname,test_params) );
  test_params.print();

  std::string spa_data_file = test_params.get<std::string>("SPA Data File");
  Int ncols    = test_params.get<Int>("ncols");
  Int nlevs    = test_params.get<Int>("nlevs");
  Int nswbands = test_params.get<Int>("nswbands");
  Int nlwbands = test_params.get<Int>("nlwbands");

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
  Int test_total_ncols = 0;
  spa_comm.all_reduce(&my_ncols,&test_total_ncols,1,MPI_SUM);
  REQUIRE(test_total_ncols == ncols);

  // Set up the set of SPA structures needed to run the test
  SPAFunc::SPAHorizInterp spa_horiz_interp;
  spa_horiz_interp.m_comm = spa_comm;
  SPAFunc::set_remap_weights_one_to_one(ncols,min_dof,dofs_gids,spa_horiz_interp);
  SPAFunc::SPATimeState     spa_time_state;
  SPAFunc::SPAData          spa_data_beg(dofs_gids.size(), nlevs+2, nswbands, nlwbands);
  SPAFunc::SPAData          spa_data_end(dofs_gids.size(), nlevs+2, nswbands, nlwbands);
  SPAFunc::SPAOutput        spa_data_out(dofs_gids.size(), nlevs,   nswbands, nlwbands);
  auto pmid_tgt = view_2d<Spack>("",dofs_gids.size(),nlevs);

  // Verify that the interpolated values match the data, since no temporal or vertical interpolation
  // should be done at this point.
  
  // Create local host views of all relevant data:
  auto hyam_h         = Kokkos::create_mirror_view(spa_data_beg.hyam);
  auto hybm_h         = Kokkos::create_mirror_view(spa_data_beg.hybm);
  // Beg data for time interpolation
  auto ps_beg         = Kokkos::create_mirror_view(spa_data_beg.PS);
  auto ccn3_beg       = Kokkos::create_mirror_view(spa_data_beg.CCN3);
  auto aer_g_sw_beg   = Kokkos::create_mirror_view(spa_data_beg.AER_G_SW);
  auto aer_ssa_sw_beg = Kokkos::create_mirror_view(spa_data_beg.AER_SSA_SW);
  auto aer_tau_sw_beg = Kokkos::create_mirror_view(spa_data_beg.AER_TAU_SW);
  auto aer_tau_lw_beg = Kokkos::create_mirror_view(spa_data_beg.AER_TAU_LW);
  // End data for time interpolation
  auto ps_end         = Kokkos::create_mirror_view(spa_data_end.PS);
  auto ccn3_end       = Kokkos::create_mirror_view(spa_data_end.CCN3);
  auto aer_g_sw_end   = Kokkos::create_mirror_view(spa_data_end.AER_G_SW);
  auto aer_ssa_sw_end = Kokkos::create_mirror_view(spa_data_end.AER_SSA_SW);
  auto aer_tau_sw_end = Kokkos::create_mirror_view(spa_data_end.AER_TAU_SW);
  auto aer_tau_lw_end = Kokkos::create_mirror_view(spa_data_end.AER_TAU_LW);
  // Output
  auto ccn3_h       = Kokkos::create_mirror_view(spa_data_out.CCN3);
  auto aer_g_sw_h   = Kokkos::create_mirror_view(spa_data_out.AER_G_SW);
  auto aer_ssa_sw_h = Kokkos::create_mirror_view(spa_data_out.AER_SSA_SW);
  auto aer_tau_sw_h = Kokkos::create_mirror_view(spa_data_out.AER_TAU_SW);
  auto aer_tau_lw_h = Kokkos::create_mirror_view(spa_data_out.AER_TAU_LW);
  
  // First initialize the start and end month data:  Set for January
  util::TimeStamp ts(1900,1,1,0,0,0);
  SPAFunc::update_spa_timestate(spa_data_file,nswbands,nlwbands,ts,spa_horiz_interp,spa_time_state,spa_data_beg,spa_data_end);

  Kokkos::deep_copy(hyam_h        ,spa_data_beg.hyam);
  Kokkos::deep_copy(hybm_h        ,spa_data_beg.hybm);
  Kokkos::deep_copy(ps_beg        ,spa_data_beg.PS);
  Kokkos::deep_copy(ccn3_beg      ,spa_data_beg.CCN3);
  Kokkos::deep_copy(aer_g_sw_beg  ,spa_data_beg.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_beg,spa_data_beg.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_beg,spa_data_beg.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_beg,spa_data_beg.AER_TAU_LW);
  Kokkos::deep_copy(ps_end        ,spa_data_end.PS);
  Kokkos::deep_copy(ccn3_end      ,spa_data_end.CCN3);
  Kokkos::deep_copy(aer_g_sw_end  ,spa_data_end.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_end,spa_data_end.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_end,spa_data_end.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_end,spa_data_end.AER_TAU_LW);

 /* ====================================================================
  * Test that spa_main is the identity when no time interpolation and
  * matching pressure profiles.
  *
  * In this section we set the target pressure profile to match the
  * expected source pressure profile that will be calculate in spa_main.
  * We test when the time stamp is both at the beginning of the month and
  * the end of the month, which should be the identity for data_out matching
  * data_beg and data_end, respectively.
  * ====================================================================*/
  // Create the pressure state.  Note, we need to create the pmid values for the actual data.  We will build it based on the PS and hya/bm
  // coordinates in the beginning data.
  auto dofs_gids_h = Kokkos::create_mirror_view(dofs_gids);
  Kokkos::deep_copy(dofs_gids_h,dofs_gids);

  auto pmid_tgt_h = Kokkos::create_mirror_view(pmid_tgt);

  // Note, hyam and hybm are padded
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack = kk / Spack::n;
      int kidx  = kk % Spack::n;
      int kpack_pad = (kk+1) / Spack::n;
      int kidx_pad  = (kk+1) % Spack::n;
      pmid_tgt_h(dof_i,kpack)[kidx] = ps_beg(dof_i)*hybm_h(kpack_pad)[kidx_pad] + P0*hyam_h(kpack_pad)[kidx_pad];
    }
  }
  int kpack_pad = (nlevs+1) / Spack::n;
  int kidx_pad  = (nlevs+1) % Spack::n;
  Kokkos::deep_copy(pmid_tgt,pmid_tgt_h);

  // Run SPA main
  SPAFunc::spa_main(spa_time_state,pmid_tgt,spa_data_beg,spa_data_end,spa_data_out,dofs_gids.size(),nlevs,nswbands,nlwbands);

  Kokkos::deep_copy(ccn3_h      , spa_data_out.CCN3);
  Kokkos::deep_copy(aer_g_sw_h  , spa_data_out.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_h, spa_data_out.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_h, spa_data_out.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_h, spa_data_out.AER_TAU_LW);

  // The output data should match the input data exactly since there is no time interpolation
  // and the pmid profile should match the constructed profile within spa_main.
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack = kk / Spack::n;
      int kidx  = kk % Spack::n;
      int kpack_pad = (kk+1) / Spack::n;
      int kidx_pad  = (kk+1) % Spack::n;
      REQUIRE(ccn3_h(dof_i,kpack)[kidx] == ccn3_beg(dof_i,kpack_pad)[kidx_pad] );
      for (int n=0;n<nswbands;n++) {
        REQUIRE(aer_g_sw_h(dof_i,n,kpack)[kidx]   == aer_g_sw_beg(dof_i,n,kpack_pad)[kidx_pad]);
        REQUIRE(aer_ssa_sw_h(dof_i,n,kpack)[kidx] == aer_ssa_sw_beg(dof_i,n,kpack_pad)[kidx_pad] );
        REQUIRE(aer_tau_sw_h(dof_i,n,kpack)[kidx] == aer_tau_sw_beg(dof_i,n,kpack_pad)[kidx_pad] );
      }
      for (int n=0;n<nlwbands;n++) {
        REQUIRE(aer_tau_lw_h(dof_i,n,kpack)[kidx] == aer_tau_lw_beg(dof_i,n,kpack_pad)[kidx_pad] );
      }
    }
  }

  /* Here we test the end of the month */
  // Add a month and recalculate.  Should now match the end of the month data.
  ts += util::days_in_month(ts.get_year(),ts.get_month())*86400;
  spa_time_state.t_now = ts.frac_of_year_in_days();

  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack = kk / Spack::n;
      int kidx  = kk % Spack::n;
      int kpack_pad = (kk+1) / Spack::n;
      int kidx_pad  = (kk+1) % Spack::n;
      pmid_tgt_h(dof_i,kpack)[kidx] = ps_end(dof_i)*hybm_h(kpack_pad)[kidx_pad] + P0*hyam_h(kpack_pad)[kidx_pad];
    }
  }
  Kokkos::deep_copy(pmid_tgt,pmid_tgt_h);
  
  SPAFunc::spa_main(spa_time_state,pmid_tgt,spa_data_beg,spa_data_end,spa_data_out,dofs_gids.size(),nlevs,nswbands,nlwbands);
  Kokkos::deep_copy(ccn3_h      , spa_data_out.CCN3);
  Kokkos::deep_copy(aer_g_sw_h  , spa_data_out.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_h, spa_data_out.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_h, spa_data_out.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_h, spa_data_out.AER_TAU_LW);

  // The output data should match the input data exactly since there is no time interpolation
  // and the pmid profile should match the constructed profile within spa_main.
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack     = kk / Spack::n;
      int kidx      = kk % Spack::n;
      int kpack_pad = (kk+1) / Spack::n;
      int kidx_pad  = (kk+1) % Spack::n;
      REQUIRE(ccn3_h(dof_i,kpack)[kidx] == ccn3_end(dof_i,kpack_pad)[kidx_pad] );
      for (int n=0;n<nswbands;n++) {
        REQUIRE(aer_g_sw_h(dof_i,n,kpack)[kidx]   == aer_g_sw_end(dof_i,n,kpack_pad)[kidx_pad]   );
        REQUIRE(aer_ssa_sw_h(dof_i,n,kpack)[kidx] == aer_ssa_sw_end(dof_i,n,kpack_pad)[kidx_pad] );
        REQUIRE(aer_tau_sw_h(dof_i,n,kpack)[kidx] == aer_tau_sw_end(dof_i,n,kpack_pad)[kidx_pad] );
      }
      for (int n=0;n<nlwbands;n++) {
        REQUIRE(aer_tau_lw_h(dof_i,n,kpack)[kidx] == aer_tau_lw_end(dof_i,n,kpack_pad)[kidx_pad] );
      }
    }
  }
 /* ====================================================================
  * Test that the time interpolation is bounded by the data:
  *
  * In the above tests we established that when the pressure profile for
  * source and data match and there is no time interpolation the output
  * data and input data match.
  *
  * Here we introduce time interpolation by picking the middle of the month
  * between beginning and ending months of data.  We set the target pressure profile
  * to match the profile expected for the source data in spa_main so that
  * only time interpolation is tested and check that all output values are
  * bounded by the min/max of the beggining/ending data.
  * ====================================================================*/
  // Add a few days and update spa data.  Make sure that the output values are not outside of the bounds
  // of the actual SPA data 
  ts += (int)(util::days_in_month(ts.get_year(),ts.get_month())*0.5)*86400;
  SPAFunc::update_spa_timestate(spa_data_file,nswbands,nlwbands,ts,spa_horiz_interp,spa_time_state,spa_data_beg,spa_data_end);
  Kokkos::deep_copy(ps_beg        ,spa_data_beg.PS);
  Kokkos::deep_copy(ccn3_beg      ,spa_data_beg.CCN3);
  Kokkos::deep_copy(aer_g_sw_beg  ,spa_data_beg.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_beg,spa_data_beg.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_beg,spa_data_beg.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_beg,spa_data_beg.AER_TAU_LW);
  Kokkos::deep_copy(ps_end        ,spa_data_end.PS);
  Kokkos::deep_copy(ccn3_end      ,spa_data_end.CCN3);
  Kokkos::deep_copy(aer_g_sw_end  ,spa_data_end.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_end,spa_data_end.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_end,spa_data_end.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_end,spa_data_end.AER_TAU_LW);
  // Create a target pressure profile to interpolate onto that matches what spa_main will compute for the source data 
  auto& t_now = spa_time_state.t_now;
  auto& t_beg = spa_time_state.t_beg_month;
  auto& t_len = spa_time_state.days_this_month;
  auto t_norm = (t_now-t_beg)/t_len;
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack     = kk / Spack::n;
      int kidx      = kk % Spack::n;
      int kpack_pad = (kk+1) / Spack::n;
      int kidx_pad  = (kk+1) % Spack::n;
      Real ps = (1.0 - t_norm) * ps_beg(dof_i) + t_norm * ps_end(dof_i);
      pmid_tgt_h(dof_i,kpack)[kidx] = ps*hybm_h(kpack_pad)[kidx_pad] + P0*hyam_h(kpack_pad)[kidx_pad];
    }
  }
  Kokkos::deep_copy(pmid_tgt,pmid_tgt_h);

  SPAFunc::spa_main(spa_time_state,pmid_tgt,spa_data_beg,spa_data_end,spa_data_out,dofs_gids.size(),nlevs,nswbands,nlwbands);
  Kokkos::deep_copy(ccn3_h      , spa_data_out.CCN3);
  Kokkos::deep_copy(aer_g_sw_h  , spa_data_out.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_h, spa_data_out.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_h, spa_data_out.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_h, spa_data_out.AER_TAU_LW);

  // Make sure the output data is within the same bounds as the input data.  We check bounds in each column and level.
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack = kk / Spack::n;
      int kidx  = kk % Spack::n;
      int kpack_pad = (kk+1) / Spack::n;
      int kidx_pad  = (kk+1) % Spack::n;
      REQUIRE(ccn3_h(dof_i,kpack)[kidx]>=std::min(ccn3_beg(dof_i,kpack_pad)[kidx_pad],ccn3_end(dof_i,kpack_pad)[kidx_pad]));
      REQUIRE(ccn3_h(dof_i,kpack)[kidx]<=std::max(ccn3_beg(dof_i,kpack_pad)[kidx_pad],ccn3_end(dof_i,kpack_pad)[kidx_pad]));
      for (int n=0;n<nswbands;n++) {
        REQUIRE(aer_g_sw_h(dof_i,n,kpack)[kidx]   >= std::min(aer_g_sw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_g_sw_end(dof_i,n,kpack_pad)[kidx_pad]) );
        REQUIRE(aer_g_sw_h(dof_i,n,kpack)[kidx]   <= std::max(aer_g_sw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_g_sw_end(dof_i,n,kpack_pad)[kidx_pad]) );

        REQUIRE(aer_ssa_sw_h(dof_i,n,kpack)[kidx] >= std::min(aer_ssa_sw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_ssa_sw_end(dof_i,n,kpack_pad)[kidx_pad]) );
        REQUIRE(aer_ssa_sw_h(dof_i,n,kpack)[kidx] <= std::max(aer_ssa_sw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_ssa_sw_end(dof_i,n,kpack_pad)[kidx_pad]) );

        REQUIRE(aer_tau_sw_h(dof_i,n,kpack)[kidx] >= std::min(aer_tau_sw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_tau_sw_end(dof_i,n,kpack_pad)[kidx_pad]) );
        REQUIRE(aer_tau_sw_h(dof_i,n,kpack)[kidx] <= std::max(aer_tau_sw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_tau_sw_end(dof_i,n,kpack_pad)[kidx_pad]) );
      }
      for (int n=0;n<nlwbands;n++) {
        REQUIRE(aer_tau_lw_h(dof_i,n,kpack)[kidx] >= std::min(aer_tau_lw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_tau_lw_end(dof_i,n,kpack_pad)[kidx_pad]) );
        REQUIRE(aer_tau_lw_h(dof_i,n,kpack)[kidx] <= std::max(aer_tau_lw_beg(dof_i,n,kpack_pad)[kidx_pad],aer_tau_lw_end(dof_i,n,kpack_pad)[kidx_pad]) );
      }
    }
  }

 /* ====================================================================
  * Test that the vertical interpolation is witin the bounds of the source data
  * on a per column basis.
  *
  * For simplicity we set data_beg = data_end which will ensure that no the
  * time interpolation returns the identity.
  *
  * We set the target pressure level to be outside of the source pressure level
  * that we expect spa_main to calculate.  We do this by bumping the surface pressure
  * for the target pressure up by 5% and halve the top-of-model pressure.
  *
  * This will require the vertical interpolation to apply the boundary conditions, which
  * should be still bounded by the max/min of the total column.
  * ====================================================================*/
  // Create a target pressure profile to interpolate onto that has a slightly higher surface pressure that the bounds.
  // This will force extrapolation.
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack     = kk / Spack::n;
      int kidx      = kk % Spack::n;
      int kpack_pad = (kk+1) / Spack::n;
      int kidx_pad  = (kk+1) % Spack::n;
      // If kk=0 then we are top-of-model, so we decrease the pressure.  Otherwise go with bigger pressure.
      Real ps   = kk==0 ? ps_beg(dof_i) : 1.05*ps_beg(dof_i);
      Real mult = kk==0 ? 0.5 : 1.0;
      pmid_tgt_h(dof_i,kpack)[kidx] = mult*(ps*hybm_h(kpack_pad)[kidx_pad] + P0*hyam_h(kpack_pad)[kidx_pad]);
    }
  }
  Kokkos::deep_copy(pmid_tgt,pmid_tgt_h);

  // Note, here we pass spa_data_beg twice, as both bounds of the input data.
  SPAFunc::spa_main(spa_time_state,pmid_tgt,spa_data_beg,spa_data_beg,spa_data_out,dofs_gids.size(),nlevs,nswbands,nlwbands);
  Kokkos::deep_copy(ccn3_h      , spa_data_out.CCN3);
  Kokkos::deep_copy(aer_g_sw_h  , spa_data_out.AER_G_SW);
  Kokkos::deep_copy(aer_ssa_sw_h, spa_data_out.AER_SSA_SW);
  Kokkos::deep_copy(aer_tau_sw_h, spa_data_out.AER_TAU_SW);
  Kokkos::deep_copy(aer_tau_lw_h, spa_data_out.AER_TAU_LW);

  // Calculate the min and max values for all spa input data for all columns,
  // note we need to compute only the interior values and ignore the padding. 
  auto ccn3_bnds       = view_2d<Real>("",2,ncols);
  auto aer_sw_g_bnds   = view_3d<Real>("",2,ncols,nswbands);
  auto aer_sw_ssa_bnds = view_3d<Real>("",2,ncols,nswbands);
  auto aer_sw_tau_bnds = view_3d<Real>("",2,ncols,nswbands);
  auto aer_lw_tau_bnds = view_3d<Real>("",2,ncols,nlwbands);
  auto ccn3_bnds_h       = Kokkos::create_mirror_view(ccn3_bnds      );
  auto aer_sw_g_bnds_h   = Kokkos::create_mirror_view(aer_sw_g_bnds  );
  auto aer_sw_ssa_bnds_h = Kokkos::create_mirror_view(aer_sw_ssa_bnds);
  auto aer_sw_tau_bnds_h = Kokkos::create_mirror_view(aer_sw_tau_bnds);
  auto aer_lw_tau_bnds_h = Kokkos::create_mirror_view(aer_lw_tau_bnds);
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    const auto& ccn3_in = ekat::subview(ccn3_beg,dof_i);
    compute_max_min(ccn3_in,1,spa_data_beg.nlevs-1,ccn3_bnds_h(0,dof_i),ccn3_bnds_h(1,dof_i));
    for (int n=0;n<nswbands;n++) {
      const auto& sw_g_in = ekat::subview(aer_g_sw_beg,dof_i,n);
      compute_max_min(sw_g_in,1,spa_data_beg.nlevs-1,aer_sw_g_bnds_h(0,dof_i,n),aer_sw_g_bnds_h(1,dof_i,n));
      const auto& sw_ssa_in = ekat::subview(aer_ssa_sw_beg,dof_i,n);
      compute_max_min(sw_ssa_in,1,spa_data_beg.nlevs-1,aer_sw_ssa_bnds_h(0,dof_i,n),aer_sw_ssa_bnds_h(1,dof_i,n));
      const auto& sw_tau_in = ekat::subview(aer_tau_sw_beg,dof_i,n);
      compute_max_min(sw_tau_in,1,spa_data_beg.nlevs-1,aer_sw_tau_bnds_h(0,dof_i,n),aer_sw_tau_bnds_h(1,dof_i,n));
    }
    for (int n=0;n<nlwbands;n++) {
      const auto& lw_tau_in = ekat::subview(aer_tau_lw_beg,dof_i,n);
      compute_max_min(lw_tau_in,1,spa_data_beg.nlevs-1,aer_lw_tau_bnds_h(0,dof_i,n),aer_lw_tau_bnds_h(1,dof_i,n));
    }
  }

  // Make sure the output data is within the same bounds as the input data.
  for (size_t dof_i=0;dof_i<dofs_gids_h.size();dof_i++) {
    for (int kk=0;kk<nlevs;kk++) {
      int kpack = kk / Spack::n;
      int kidx  = kk % Spack::n;
      REQUIRE(ccn3_h(dof_i,kpack)[kidx]>=ccn3_bnds_h(0,dof_i));
      REQUIRE(ccn3_h(dof_i,kpack)[kidx]<=ccn3_bnds_h(1,dof_i));
      for (int n=0;n<nswbands;n++) {
        REQUIRE(aer_g_sw_h(dof_i,n,kpack)[kidx]   >= aer_sw_g_bnds_h  (0,dof_i,n) );
        REQUIRE(aer_ssa_sw_h(dof_i,n,kpack)[kidx] >= aer_sw_ssa_bnds_h(0,dof_i,n) );
        REQUIRE(aer_tau_sw_h(dof_i,n,kpack)[kidx] >= aer_sw_tau_bnds_h(0,dof_i,n) );
        REQUIRE(aer_g_sw_h(dof_i,n,kpack)[kidx]   <= aer_sw_g_bnds_h  (1,dof_i,n) );
        REQUIRE(aer_ssa_sw_h(dof_i,n,kpack)[kidx] <= aer_sw_ssa_bnds_h(1,dof_i,n) );
        REQUIRE(aer_tau_sw_h(dof_i,n,kpack)[kidx] <= aer_sw_tau_bnds_h(1,dof_i,n) );
      }
      for (int n=0;n<nlwbands;n++) {
        REQUIRE(aer_tau_lw_h(dof_i,n,kpack)[kidx] >= aer_lw_tau_bnds_h(0,dof_i,n) );
        REQUIRE(aer_tau_lw_h(dof_i,n,kpack)[kidx] <= aer_lw_tau_bnds_h(1,dof_i,n) );
      }
    }
  }

 /* ====================================================================
  *                         DONE WITH TESTS
  * ====================================================================*/
  // All Done 
  scorpio::eam_pio_finalize();
} // run_property

// Helper function
void compute_max_min(const view_1d<const Spack>& input, const int start, const int end, Real& min, Real& max)
{
  int kpack, kidx;
  // Initialize min and max with first entry in input
  kpack = start / Spack::n;
  kidx  = start % Spack::n;
  min = input(kpack)[kidx];
  max = input(kpack)[kidx];
  // Now go through the rest of the view
  for (int kk = start+1; kk<=end; kk++) {
    kpack = kk / Spack::n;
    kidx  = kk % Spack::n;
    min = std::min(min,input(kpack)[kidx]);
    max = std::max(max,input(kpack)[kidx]);
  }
}

} // namespace

