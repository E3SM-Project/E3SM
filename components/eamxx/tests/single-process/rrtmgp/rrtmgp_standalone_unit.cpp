#include <catch2/catch.hpp>

// The atmosphere driver
#include "control/atmosphere_driver.hpp"

// Other rrtmgp specific code needed specifically for this test
#include "physics/rrtmgp/rrtmgp_test_utils.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/eamxx_rrtmgp_process_interface.hpp"
#include "physics/register_physics.hpp"
#include "physics/share/physics_constants.hpp"

// scream share headers
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/scream_common_physics_functions.hpp"

// EKAT headers
#include <ekat/ekat_parse_yaml_file.hpp>
#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/util/ekat_test_utils.hpp>

// RRTMGP and YAKL
#include <mo_gas_concentrations.h>
#include <mo_garand_atmos_io.h>
#ifdef RRTMGP_ENABLE_YAKL
#include <YAKL.h>
#endif

// System headers
#include <iostream>
#include <cmath>

/*
 * Run standalone test problem for RRTMGP and compare with baseline
 */

namespace {

using namespace scream;
using namespace scream::control;
using KT = KokkosTypes<DefaultDevice>;
using ExeSpace = KT::ExeSpace;
using MemberType = KT::MemberType;
using PC = scream::physics::Constants<Real>;

/*
 * Run standalone test through SCREAM driver this time
 */
#ifdef RRTMGP_ENABLE_YAKL
TEST_CASE("rrtmgp_scream_standalone", "") {
  // Get baseline name (needs to be passed as an arg)
  std::string inputfile = ekat::TestSession::get().params.at("inputfile");
  std::string baseline = ekat::TestSession::get().params.at("baseline");

  // Check if files exists
  REQUIRE(rrtmgpTest::file_exists(inputfile.c_str()));
  REQUIRE(rrtmgpTest::file_exists(baseline.c_str()));

  // Initialize yakl
  if(!yakl::isInitialized()) { yakl::init(); }

  // Read reference fluxes from baseline file
  real2d sw_flux_up_ref;
  real2d sw_flux_dn_ref;
  real2d sw_flux_dn_dir_ref;
  real2d lw_flux_up_ref;
  real2d lw_flux_dn_ref;
  rrtmgpTest::read_fluxes(baseline, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

  // Load ad parameter list
  std::string fname = "input_unit.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);
  // Create a MPI communicator
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  register_physics ();
  register_mesh_free_grids_manager();

  // Create the driver
  AtmosphereDriver ad;

  // Dummy timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});

  // Initialize the driver
  ad.initialize(atm_comm, ad_params, time);

  /*
   * Setup the dummy problem and overwrite default initial conditions
   */

  // Get dimension sizes from the field manager
  const auto& grid = ad.get_grids_manager()->get_grid("Point Grid");
  const auto& field_mgr = *ad.get_field_mgr(grid->name());
  int ncol  = grid->get_num_local_dofs();
  int nlay  = grid->get_num_vertical_levels();

  // In this test, we need to hack lat/lon. But the fields we get
  // from the grid are read-only. Therefore, hack a bit, and cast
  // away constness. It's bad, but it's only for this unit test
  auto clat = grid->get_geometry_data("lat").get_view<const Real*>();
  auto clon = grid->get_geometry_data("lon").get_view<const Real*>();
  auto lat = const_cast<Real*>(clat.data());
  auto lon = const_cast<Real*>(clon.data());

  // Get number of shortwave bands and number of gases from RRTMGP
  int ngas     =   8;  // TODO: get this intelligently

  // Make sure we have the right dimension sizes
  REQUIRE(nlay == static_cast<int>(sw_flux_up_ref.dimension[1])-1);

  // Create yakl arrays to store the input data
  auto p_lay = real2d("p_lay", ncol, nlay);
  auto t_lay = real2d("t_lay", ncol, nlay);
  auto p_del = real2d("p_del", ncol, nlay);
  auto p_lev = real2d("p_lev", ncol, nlay+1);
  auto t_lev = real2d("t_lev", ncol, nlay+1);
  auto sfc_alb_dir_vis = real1d("sfc_alb_dir_vis", ncol);
  auto sfc_alb_dir_nir = real1d("sfc_alb_dir_nir", ncol);
  auto sfc_alb_dif_vis = real1d("sfc_alb_dif_vis", ncol);
  auto sfc_alb_dif_nir = real1d("sfc_alb_dif_nir", ncol);
  auto lwp = real2d("lwp", ncol, nlay);
  auto iwp = real2d("iwp", ncol, nlay);
  auto rel = real2d("rel", ncol, nlay);
  auto rei = real2d("rei", ncol, nlay);
  auto cld = real2d("cld", ncol, nlay);
  auto mu0 = real1d("mu0", ncol);
  auto gas_vmr = real3d("gas_vmr", ncol, nlay, ngas);

  // Setup dummy problem; need to use tmp arrays with ncol_all size
  auto ncol_all = grid->get_num_global_dofs();
  auto p_lay_all = real2d("p_lay", ncol_all, nlay);
  auto t_lay_all = real2d("t_lay", ncol_all, nlay);
  auto p_lev_all = real2d("p_lev", ncol_all, nlay+1);
  auto t_lev_all = real2d("t_lev", ncol_all, nlay+1);
  auto sfc_alb_dir_vis_all = real1d("sfc_alb_dir_vis", ncol_all);
  auto sfc_alb_dir_nir_all = real1d("sfc_alb_dir_nir", ncol_all);
  auto sfc_alb_dif_vis_all = real1d("sfc_alb_dif_vis", ncol_all);
  auto sfc_alb_dif_nir_all = real1d("sfc_alb_dif_nir", ncol_all);
  auto lwp_all = real2d("lwp", ncol_all, nlay);
  auto iwp_all = real2d("iwp", ncol_all, nlay);
  auto rel_all = real2d("rel", ncol_all, nlay);
  auto rei_all = real2d("rei", ncol_all, nlay);
  auto cld_all = real2d("cld", ncol_all, nlay);
  auto mu0_all = real1d("mu0", ncol_all);
  // Read in dummy Garand atmosphere; if this were an actual model simulation,
  // these would be passed as inputs to the driver
  // NOTE: set ncol to size of col_flx dimension in the input file. This is so
  // that we can compare to the reference data provided in that file. Note that
  // this will copy the first column of the input data (the first profile) ncol
  // times. We will then fill some fraction of these columns with clouds for
  // the test problem.
  GasConcs gas_concs;
  real2d col_dry;
  read_atmos(inputfile, p_lay_all, t_lay_all, p_lev_all, t_lev_all, gas_concs, col_dry, ncol_all);
  // Setup dummy problem; need to use tmp arrays with ncol_all size
  rrtmgpTest::dummy_atmos(
    inputfile, ncol_all, p_lay_all, t_lay_all,
    sfc_alb_dir_vis_all, sfc_alb_dir_nir_all,
    sfc_alb_dif_vis_all, sfc_alb_dif_nir_all,
    mu0_all,
    lwp_all, iwp_all, rel_all, rei_all, cld_all
                          );
  // Populate our local input arrays with the proper columns, based on our rank
  int irank = atm_comm.rank();
  yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<1>(ncol), YAKL_LAMBDA(int icol) {
      auto icol_all = ncol * irank + icol;
      sfc_alb_dir_vis(icol) = sfc_alb_dir_vis_all(icol_all);
      sfc_alb_dir_nir(icol) = sfc_alb_dir_nir_all(icol_all);
      sfc_alb_dif_vis(icol) = sfc_alb_dif_vis_all(icol_all);
      sfc_alb_dif_nir(icol) = sfc_alb_dif_nir_all(icol_all);
      mu0(icol) = mu0_all(icol_all);
    });
  yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
      auto icol_all = ncol * irank + icol;
      p_lay(icol,ilay) = p_lay_all(icol_all,ilay);
      t_lay(icol,ilay) = t_lay_all(icol_all,ilay);
      lwp(icol,ilay) = lwp_all(icol_all,ilay);
      iwp(icol,ilay) = iwp_all(icol_all,ilay);
      rel(icol,ilay) = rel_all(icol_all,ilay);
      rei(icol,ilay) = rei_all(icol_all,ilay);
      cld(icol,ilay) = cld_all(icol_all,ilay);
    });
  yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay+1,ncol), YAKL_LAMBDA(int ilay, int icol) {
      auto icol_all = ncol * irank + icol;
      p_lev(icol,ilay) = p_lev_all(icol_all,ilay);
      t_lev(icol,ilay) = t_lev_all(icol_all,ilay);
    });
  // Free temporary variables
  p_lay_all.deallocate();
  t_lay_all.deallocate();
  p_lev_all.deallocate();
  t_lev_all.deallocate();
  sfc_alb_dir_vis_all.deallocate();
  sfc_alb_dir_nir_all.deallocate();
  sfc_alb_dif_vis_all.deallocate();
  sfc_alb_dif_nir_all.deallocate();
  lwp_all.deallocate();
  iwp_all.deallocate();
  rel_all.deallocate();
  rei_all.deallocate();
  cld_all.deallocate();
  mu0_all.deallocate();

  // Need to calculate a dummy pseudo_density for our test problem
  yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol) {
      p_del(icol,ilay) = abs(p_lev(icol,ilay+1) - p_lev(icol,ilay));
    });

  // We do not want to pass lwp and iwp through the FM, so back out qc and qi for this test
  // NOTE: test problem provides lwp/iwp in g/m2, not kg/m2! Factor of 1e3 here!
  auto qc = real2d("qc", ncol, nlay);
  auto nc = real2d("nc", ncol, nlay);
  auto qi = real2d("qi", ncol, nlay);
  yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay, ncol), YAKL_LAMBDA(int ilay, int icol) {
      qc(icol,ilay) = 1e-3 * lwp(icol,ilay) * cld(icol,ilay) * PC::gravit / p_del(icol,ilay);
      qi(icol,ilay) = 1e-3 * iwp(icol,ilay) * cld(icol,ilay) * PC::gravit / p_del(icol,ilay);
    });

  // Copy gases from gas_concs to gas_vmr array
  yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<3>(ncol,nlay,ngas), YAKL_LAMBDA(int icol, int ilay, int igas) {
      auto icol_all = ncol * irank + icol;
      gas_vmr(icol,ilay,igas) = gas_concs.concs(icol_all,ilay,igas);
    });
  gas_concs.reset();

  // Before running, make a copy of T_mid so we can see changes
  auto T_mid0 = real2d("T_mid0", ncol, nlay);
  t_lay.deep_copy_to(T_mid0);

  // Grab views from field manager and copy in values from yakl arrays. Making
  // copies is necessary since the yakl::Array take in data arranged with ncol
  // as the fastest index, but the field manager expects the 2nd dimension as
  // the fastest index.
  auto d_pmid = field_mgr.get_field("p_mid").get_view<Real**>();
  auto d_tmid = field_mgr.get_field("T_mid").get_view<Real**>();
  auto d_pint = field_mgr.get_field("p_int").get_view<Real**>();
  auto d_pdel = field_mgr.get_field("pseudo_density").get_view<Real**>();
  auto d_sfc_alb_dir_vis = field_mgr.get_field("sfc_alb_dir_vis").get_view<Real*>();
  auto d_sfc_alb_dir_nir = field_mgr.get_field("sfc_alb_dir_nir").get_view<Real*>();
  auto d_sfc_alb_dif_vis = field_mgr.get_field("sfc_alb_dif_vis").get_view<Real*>();
  auto d_sfc_alb_dif_nir = field_mgr.get_field("sfc_alb_dif_nir").get_view<Real*>();
  auto d_surf_lw_flux_up = field_mgr.get_field("surf_lw_flux_up").get_view<Real*>();
  auto d_qc = field_mgr.get_field("qc").get_view<Real**>();
  auto d_nc = field_mgr.get_field("nc").get_view<Real**>();
  auto d_qi = field_mgr.get_field("qi").get_view<Real**>();
  auto d_rel = field_mgr.get_field("eff_radius_qc").get_view<Real**>();
  auto d_rei = field_mgr.get_field("eff_radius_qi").get_view<Real**>();
  auto d_cld = field_mgr.get_field("cldfrac_tot").get_view<Real**>();

  auto d_qv  = field_mgr.get_field("qv").get_view<Real**>();
  auto d_co2 = field_mgr.get_field("co2_volume_mix_ratio").get_view<Real**>();
  auto d_o3  = field_mgr.get_field("o3_volume_mix_ratio").get_view<Real**>();
  auto d_n2o = field_mgr.get_field("n2o_volume_mix_ratio").get_view<Real**>();
  auto d_co  = field_mgr.get_field("co_volume_mix_ratio").get_view<Real**>();
  auto d_ch4 = field_mgr.get_field("ch4_volume_mix_ratio").get_view<Real**>();
  auto d_o2  = field_mgr.get_field("o2_volume_mix_ratio").get_view<Real**>();
  auto d_n2  = field_mgr.get_field("n2_volume_mix_ratio").get_view<Real**>();

  // Gather molecular weights of all the active gases in the test for conversion
  // to mass-mixing-ratio.
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();

        // Set lat and lon to single value for just this test:
        // Note, these values will ensure that the cosine zenith
        // angle will end up matching the constant value meant for
        // the test, which is 0.86
        lat[i] = 5.224000000000;
        lon[i] = 167.282000000000;

        d_sfc_alb_dir_vis(i) = sfc_alb_dir_vis(i+1);
        d_sfc_alb_dir_nir(i) = sfc_alb_dir_nir(i+1);
        d_sfc_alb_dif_vis(i) = sfc_alb_dif_vis(i+1);
        d_sfc_alb_dif_nir(i) = sfc_alb_dif_nir(i+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
            d_pmid(i,k) = p_lay(i+1,k+1);
            d_tmid(i,k) = t_lay(i+1,k+1);
            d_pdel(i,k) = p_del(i+1,k+1);
            d_qc(i,k)   = qc(i+1,k+1);
            d_nc(i,k)   = nc(i+1,k+1);
            d_qi(i,k)   = qi(i+1,k+1);
            d_rel(i,k)  = rel(i+1,k+1);
            d_rei(i,k)  = rei(i+1,k+1);
            d_cld(i,k)  = cld(i+1,k+1);
            d_pint(i,k) = p_lev(i+1,k+1);
            // qv specified as a wet mixing ratio
            Real qv_dry = gas_vmr(i+1,k+1,1)*PC::ep_2;
            Real qv_wet = qv_dry/(1.0+qv_dry);
            d_qv(i,k)  = qv_wet;
            // rest of active gases are specified as volume mixing ratios
            d_co2(i,k) = gas_vmr(i+1,k+1,2);
            d_o3(i,k)  = gas_vmr(i+1,k+1,3);
            d_n2o(i,k) = gas_vmr(i+1,k+1,4);
            d_co(i,k)  = gas_vmr(i+1,k+1,5);
            d_ch4(i,k) = gas_vmr(i+1,k+1,6);
            d_o2(i,k)  = gas_vmr(i+1,k+1,7);
            d_n2(i,k)  = gas_vmr(i+1,k+1,8);
          });

        d_pint(i,nlay) = p_lev(i+1,nlay+1);

        // Compute surface flux from surface temperature
        auto ibot = (p_lay(1,1) > p_lay(1,nlay)) ? 1 : nlay+1;
        d_surf_lw_flux_up(i) = PC::stebol * pow(t_lev(i+1,ibot), 4.0);
      });
  }
  Kokkos::fence();

  // Run driver
  ad.run(300);

  // Check values; The correct values have been stored in the field manager, we need to
  // copy back to YAKL::Array.
  auto d_sw_flux_up = field_mgr.get_field("sw_flux_up").get_view<Real**>();
  auto d_sw_flux_dn = field_mgr.get_field("sw_flux_dn").get_view<Real**>();
  auto d_sw_flux_dn_dir = field_mgr.get_field("sw_flux_dn_dir").get_view<Real**>();
  auto d_lw_flux_up = field_mgr.get_field("lw_flux_up").get_view<Real**>();
  auto d_lw_flux_dn = field_mgr.get_field("lw_flux_dn").get_view<Real**>();
  auto sw_flux_up_test = real2d("sw_flux_up_test", ncol, nlay+1);
  auto sw_flux_dn_test = real2d("sw_flux_dn_test", ncol, nlay+1);
  auto sw_flux_dn_dir_test = real2d("sw_flux_dn_dir_test",  ncol, nlay+1);
  auto lw_flux_up_test = real2d("lw_flux_up_test", ncol, nlay+1);
  auto lw_flux_dn_test = real2d("lw_flux_dn_test", ncol, nlay+1);
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
        const int i = team.league_rank();

        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay+1), [&] (const int& k) {
            if (k < nlay) t_lay(i+1,k+1) = d_tmid(i,k);

            sw_flux_up_test(i+1,k+1)     = d_sw_flux_up(i,k);
            sw_flux_dn_test(i+1,k+1)     = d_sw_flux_dn(i,k);
            sw_flux_dn_dir_test(i+1,k+1) = d_sw_flux_dn_dir(i,k);
            lw_flux_up_test(i+1,k+1)     = d_lw_flux_up(i,k);
            lw_flux_dn_test(i+1,k+1)     = d_lw_flux_dn(i,k);
          });
      });
  }
  Kokkos::fence();

  // Dumb check to verify that we did indeed update temperature
  REQUIRE(t_lay.createHostCopy()(1,1) != T_mid0.createHostCopy()(1,1));
  T_mid0.deallocate();

  // Make sure fluxes from field manager that were calculated in AD call of RRTMGP match reference fluxes from input file
  // We use all_close here instead of all_equals because we are only able
  // to approximate the solar zenith angle used in the RRTMGP clear-sky
  // test problem with our trial and error lat/lon values, so fluxes will
  // be slightly off. We just verify that they are all "close" here, within
  // some tolerance. Computation of level temperatures from midpoints is
  // also unable to exactly reproduce the values in the test problem, so
  // computed fluxes will be further off from the reference calculation.
  //
  // We need to create local copies with only the columns specific to our rank in case we have split columns over multiple ranks
  auto sw_flux_up_loc     = real2d("sw_flux_up_loc"    , ncol, nlay+1);
  auto sw_flux_dn_loc     = real2d("sw_flux_dn_loc"    , ncol, nlay+1);
  auto sw_flux_dn_dir_loc = real2d("sw_flux_dn_dir_loc", ncol, nlay+1);
  auto lw_flux_up_loc     = real2d("lw_flux_up_loc"    , ncol, nlay+1);
  auto lw_flux_dn_loc     = real2d("lw_flux_dn_loc"    , ncol, nlay+1);
  yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay+1,ncol), YAKL_LAMBDA(int ilay, int icol) {
      auto icol_all = ncol * irank + icol;
      sw_flux_up_loc(icol,ilay) = sw_flux_up_ref(icol_all,ilay);
      sw_flux_dn_loc(icol,ilay) = sw_flux_dn_ref(icol_all,ilay);
      sw_flux_dn_dir_loc(icol,ilay) = sw_flux_dn_dir_ref(icol_all,ilay);
      lw_flux_up_loc(icol,ilay) = lw_flux_up_ref(icol_all,ilay);
      lw_flux_dn_loc(icol,ilay) = lw_flux_dn_ref(icol_all,ilay);
    });
  REQUIRE(rrtmgpTest::all_close(sw_flux_up_loc    , sw_flux_up_test    , 1.0));
  REQUIRE(rrtmgpTest::all_close(sw_flux_dn_loc    , sw_flux_dn_test    , 1.0));
  REQUIRE(rrtmgpTest::all_close(sw_flux_dn_dir_loc, sw_flux_dn_dir_test, 1.0));
  REQUIRE(rrtmgpTest::all_close(lw_flux_up_loc    , lw_flux_up_test    , 1.0));
  REQUIRE(rrtmgpTest::all_close(lw_flux_dn_loc    , lw_flux_dn_test    , 1.0));

  // Deallocate YAKL arrays
  p_lay.deallocate();
  t_lay.deallocate();
  p_del.deallocate();
  p_lev.deallocate();
  t_lev.deallocate();
  sfc_alb_dir_vis.deallocate();
  sfc_alb_dir_nir.deallocate();
  sfc_alb_dif_vis.deallocate();
  sfc_alb_dif_nir.deallocate();
  lwp.deallocate();
  iwp.deallocate();
  rel.deallocate();
  rei.deallocate();
  cld.deallocate();
  qc.deallocate();
  nc.deallocate();
  qi.deallocate();
  mu0.deallocate();
  gas_vmr.deallocate();
  sw_flux_up_test.deallocate();
  sw_flux_dn_test.deallocate();
  sw_flux_dn_dir_test.deallocate();
  lw_flux_up_test.deallocate();
  lw_flux_dn_test.deallocate();
  sw_flux_up_ref.deallocate();
  sw_flux_dn_ref.deallocate();
  sw_flux_dn_dir_ref.deallocate();
  lw_flux_up_ref.deallocate();
  lw_flux_dn_ref.deallocate();
  sw_flux_up_loc.deallocate();
  sw_flux_dn_loc.deallocate();
  sw_flux_dn_dir_loc.deallocate();
  lw_flux_up_loc.deallocate();
  lw_flux_dn_loc.deallocate();
  col_dry.deallocate();

  // Finalize the driver. YAKL will be finalized inside
  // RRTMGPRadiation::finalize_impl after RRTMGP has had the
  // opportunity to deallocate all it's arrays.
  ad.finalize();
}
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
TEST_CASE("rrtmgp_scream_standalone_k", "") {
#ifdef RRTMGP_LAYOUT_LEFT
  using layout_t   = Kokkos::LayoutLeft;
#else
  using layout_t   = typename ekat::KokkosTypes<DefaultDevice>::Layout;
#endif
  using interface_t = scream::rrtmgp::rrtmgp_interface<Real, layout_t, DefaultDevice>;
  using MDRP = interface_t::MDRP;
  using utils_t = rrtmgpTest::rrtmgp_test_utils<Real, layout_t, DefaultDevice>;
  using real1dk = interface_t::view_t<Real*>;
  using real2dk = interface_t::view_t<Real**>;
  using real3dk = interface_t::view_t<Real***>;

  // Get baseline name (needs to be passed as an arg)
  std::string inputfile = ekat::TestSession::get().params.at("rrtmgp_inputfile");
  std::string baseline = ekat::TestSession::get().params.at("rrtmgp_baseline");

  // Check if files exists
  REQUIRE(rrtmgpTest::file_exists(inputfile.c_str()));
  REQUIRE(rrtmgpTest::file_exists(baseline.c_str()));

  // Initialize yakl
  scream::init_kls();

  // Read reference fluxes from baseline file
  real2dk sw_flux_up_ref;
  real2dk sw_flux_dn_ref;
  real2dk sw_flux_dn_dir_ref;
  real2dk lw_flux_up_ref;
  real2dk lw_flux_dn_ref;
  utils_t::read_fluxes(baseline, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

  // Load ad parameter list
  std::string fname = "input_unit.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);
  // Create a MPI communicator
  ekat::Comm atm_comm (MPI_COMM_WORLD);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  register_physics ();
  register_mesh_free_grids_manager();

  // Create the driver
  AtmosphereDriver ad;

  // Dummy timestamp
  util::TimeStamp time ({2000,1,1},{0,0,0});

  // Initialize the driver
  ad.initialize(atm_comm, ad_params, time);

  /*
   * Setup the dummy problem and overwrite default initial conditions
   */

  // Get dimension sizes from the field manager
  const auto& grid = ad.get_grids_manager()->get_grid("Point Grid");
  const auto& field_mgr = *ad.get_field_mgr(grid->name());
  int ncol  = grid->get_num_local_dofs();
  int nlay  = grid->get_num_vertical_levels();

  // In this test, we need to hack lat/lon. But the fields we get
  // from the grid are read-only. Therefore, hack a bit, and cast
  // away constness. It's bad, but it's only for this unit test
  auto clat = grid->get_geometry_data("lat").get_view<const Real*>();
  auto clon = grid->get_geometry_data("lon").get_view<const Real*>();
  auto lat = const_cast<Real*>(clat.data());
  auto lon = const_cast<Real*>(clon.data());

  // Get number of shortwave bands and number of gases from RRTMGP
  int ngas     =   8;  // TODO: get this intelligently

  // Make sure we have the right dimension sizes
  REQUIRE(nlay == static_cast<int>(sw_flux_up_ref.extent(1))-1);

  // Create yakl arrays to store the input data
  auto p_lay = real2dk("p_lay", ncol, nlay);
  auto t_lay = real2dk("t_lay", ncol, nlay);
  auto p_del = real2dk("p_del", ncol, nlay);
  auto p_lev = real2dk("p_lev", ncol, nlay+1);
  auto t_lev = real2dk("t_lev", ncol, nlay+1);
  auto sfc_alb_dir_vis = real1dk("sfc_alb_dir_vis", ncol);
  auto sfc_alb_dir_nir = real1dk("sfc_alb_dir_nir", ncol);
  auto sfc_alb_dif_vis = real1dk("sfc_alb_dif_vis", ncol);
  auto sfc_alb_dif_nir = real1dk("sfc_alb_dif_nir", ncol);
  auto lwp = real2dk("lwp", ncol, nlay);
  auto iwp = real2dk("iwp", ncol, nlay);
  auto rel = real2dk("rel", ncol, nlay);
  auto rei = real2dk("rei", ncol, nlay);
  auto cld = real2dk("cld", ncol, nlay);
  auto mu0 = real1dk("mu0", ncol);
  auto gas_vmr = real3dk("gas_vmr", ncol, nlay, ngas);

  // Setup dummy problem; need to use tmp arrays with ncol_all size
  auto ncol_all = grid->get_num_global_dofs();
  auto p_lay_all = real2dk("p_lay", ncol_all, nlay);
  auto t_lay_all = real2dk("t_lay", ncol_all, nlay);
  auto p_lev_all = real2dk("p_lev", ncol_all, nlay+1);
  auto t_lev_all = real2dk("t_lev", ncol_all, nlay+1);
  auto sfc_alb_dir_vis_all = real1dk("sfc_alb_dir_vis", ncol_all);
  auto sfc_alb_dir_nir_all = real1dk("sfc_alb_dir_nir", ncol_all);
  auto sfc_alb_dif_vis_all = real1dk("sfc_alb_dif_vis", ncol_all);
  auto sfc_alb_dif_nir_all = real1dk("sfc_alb_dif_nir", ncol_all);
  auto lwp_all = real2dk("lwp", ncol_all, nlay);
  auto iwp_all = real2dk("iwp", ncol_all, nlay);
  auto rel_all = real2dk("rel", ncol_all, nlay);
  auto rei_all = real2dk("rei", ncol_all, nlay);
  auto cld_all = real2dk("cld", ncol_all, nlay);
  auto mu0_all = real1dk("mu0", ncol_all);
  // Read in dummy Garand atmosphere; if this were an actual model simulation,
  // these would be passed as inputs to the driver
  // NOTE: set ncol to size of col_flx dimension in the input file. This is so
  // that we can compare to the reference data provided in that file. Note that
  // this will copy the first column of the input data (the first profile) ncol
  // times. We will then fill some fraction of these columns with clouds for
  // the test problem.
  GasConcsK<Real, Kokkos::LayoutRight, DefaultDevice> gas_concs;
  real2dk col_dry;
  read_atmos(inputfile, p_lay_all, t_lay_all, p_lev_all, t_lev_all, gas_concs, col_dry, ncol_all);
  // Setup dummy problem; need to use tmp arrays with ncol_all size
  utils_t::dummy_atmos(
    inputfile, ncol_all, p_lay_all, t_lay_all,
    sfc_alb_dir_vis_all, sfc_alb_dir_nir_all,
    sfc_alb_dif_vis_all, sfc_alb_dif_nir_all,
    mu0_all,
    lwp_all, iwp_all, rel_all, rei_all, cld_all
                          );
  // Populate our local input arrays with the proper columns, based on our rank
  int irank = atm_comm.rank();
  Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    auto icol_all = ncol * irank + icol;
    sfc_alb_dir_vis(icol) = sfc_alb_dir_vis_all(icol_all);
    sfc_alb_dir_nir(icol) = sfc_alb_dir_nir_all(icol_all);
    sfc_alb_dif_vis(icol) = sfc_alb_dif_vis_all(icol_all);
    sfc_alb_dif_nir(icol) = sfc_alb_dif_nir_all(icol_all);
    mu0(icol) = mu0_all(icol_all);
  });
  Kokkos::parallel_for(MDRP::template get<2>({nlay,ncol}), KOKKOS_LAMBDA(int ilay, int icol) {
    auto icol_all = ncol * irank + icol;
    p_lay(icol,ilay) = p_lay_all(icol_all,ilay);
    t_lay(icol,ilay) = t_lay_all(icol_all,ilay);
    lwp(icol,ilay) = lwp_all(icol_all,ilay);
    iwp(icol,ilay) = iwp_all(icol_all,ilay);
    rel(icol,ilay) = rel_all(icol_all,ilay);
    rei(icol,ilay) = rei_all(icol_all,ilay);
    cld(icol,ilay) = cld_all(icol_all,ilay);
  });
  Kokkos::parallel_for(MDRP::template get<2>({nlay+1,ncol}), KOKKOS_LAMBDA(int ilay, int icol) {
    auto icol_all = ncol * irank + icol;
    p_lev(icol,ilay) = p_lev_all(icol_all,ilay);
    t_lev(icol,ilay) = t_lev_all(icol_all,ilay);
  });

  // Need to calculate a dummy pseudo_density for our test problem
  Kokkos::parallel_for(MDRP::template get<2>({nlay,ncol}), KOKKOS_LAMBDA(int ilay, int icol) {
    p_del(icol,ilay) = abs(p_lev(icol,ilay+1) - p_lev(icol,ilay));
  });

  // We do not want to pass lwp and iwp through the FM, so back out qc and qi for this test
  // NOTE: test problem provides lwp/iwp in g/m2, not kg/m2! Factor of 1e3 here!
  auto qc = real2dk("qc", ncol, nlay);
  auto nc = real2dk("nc", ncol, nlay);
  auto qi = real2dk("qi", ncol, nlay);
  Kokkos::parallel_for(MDRP::template get<2>({nlay,ncol}), KOKKOS_LAMBDA(int ilay, int icol) {
    qc(icol,ilay) = 1e-3 * lwp(icol,ilay) * cld(icol,ilay) * PC::gravit / p_del(icol,ilay);
    qi(icol,ilay) = 1e-3 * iwp(icol,ilay) * cld(icol,ilay) * PC::gravit / p_del(icol,ilay);
  });

  // Copy gases from gas_concs to gas_vmr array
  Kokkos::parallel_for(MDRP::template get<3>({ncol,nlay,ngas}), KOKKOS_LAMBDA(int icol, int ilay, int igas) {
    auto icol_all = ncol * irank + icol;
    gas_vmr(icol,ilay,igas) = gas_concs.concs(icol_all,ilay,igas);
  });
  gas_concs.reset();

  // Before running, make a copy of T_mid so we can see changes
  auto T_mid0 = real2dk("T_mid0", ncol, nlay);
  Kokkos::deep_copy(T_mid0, t_lay);

  // Grab views from field manager and copy in values from yakl arrays. Making
  // copies is necessary since the yakl::Array take in data arranged with ncol
  // as the fastest index, but the field manager expects the 2nd dimension as
  // the fastest index.
  auto d_pmid = field_mgr.get_field("p_mid").get_view<Real**>();
  auto d_tmid = field_mgr.get_field("T_mid").get_view<Real**>();
  auto d_pint = field_mgr.get_field("p_int").get_view<Real**>();
  auto d_pdel = field_mgr.get_field("pseudo_density").get_view<Real**>();
  auto d_sfc_alb_dir_vis = field_mgr.get_field("sfc_alb_dir_vis").get_view<Real*>();
  auto d_sfc_alb_dir_nir = field_mgr.get_field("sfc_alb_dir_nir").get_view<Real*>();
  auto d_sfc_alb_dif_vis = field_mgr.get_field("sfc_alb_dif_vis").get_view<Real*>();
  auto d_sfc_alb_dif_nir = field_mgr.get_field("sfc_alb_dif_nir").get_view<Real*>();
  auto d_surf_lw_flux_up = field_mgr.get_field("surf_lw_flux_up").get_view<Real*>();
  auto d_qc = field_mgr.get_field("qc").get_view<Real**>();
  auto d_nc = field_mgr.get_field("nc").get_view<Real**>();
  auto d_qi = field_mgr.get_field("qi").get_view<Real**>();
  auto d_rel = field_mgr.get_field("eff_radius_qc").get_view<Real**>();
  auto d_rei = field_mgr.get_field("eff_radius_qi").get_view<Real**>();
  auto d_cld = field_mgr.get_field("cldfrac_tot").get_view<Real**>();

  auto d_qv  = field_mgr.get_field("qv").get_view<Real**>();
  auto d_co2 = field_mgr.get_field("co2_volume_mix_ratio").get_view<Real**>();
  auto d_o3  = field_mgr.get_field("o3_volume_mix_ratio").get_view<Real**>();
  auto d_n2o = field_mgr.get_field("n2o_volume_mix_ratio").get_view<Real**>();
  auto d_co  = field_mgr.get_field("co_volume_mix_ratio").get_view<Real**>();
  auto d_ch4 = field_mgr.get_field("ch4_volume_mix_ratio").get_view<Real**>();
  auto d_o2  = field_mgr.get_field("o2_volume_mix_ratio").get_view<Real**>();
  auto d_n2  = field_mgr.get_field("n2_volume_mix_ratio").get_view<Real**>();

  // Gather molecular weights of all the active gases in the test for conversion
  // to mass-mixing-ratio.
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();

      // Set lat and lon to single value for just this test:
      // Note, these values will ensure that the cosine zenith
      // angle will end up matching the constant value meant for
      // the test, which is 0.86
      lat[i] = 5.224000000000;
      lon[i] = 167.282000000000;

      d_sfc_alb_dir_vis(i) = sfc_alb_dir_vis(i);
      d_sfc_alb_dir_nir(i) = sfc_alb_dir_nir(i);
      d_sfc_alb_dif_vis(i) = sfc_alb_dif_vis(i);
      d_sfc_alb_dif_nir(i) = sfc_alb_dif_nir(i);
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay), [&] (const int& k) {
        d_pmid(i,k) = p_lay(i,k);
        d_tmid(i,k) = t_lay(i,k);
        d_pdel(i,k) = p_del(i,k);
        d_qc(i,k)   = qc(i,k);
        d_nc(i,k)   = nc(i,k);
        d_qi(i,k)   = qi(i,k);
        d_rel(i,k)  = rel(i,k);
        d_rei(i,k)  = rei(i,k);
        d_cld(i,k)  = cld(i,k);
        d_pint(i,k) = p_lev(i,k);
        // qv specified as a wet mixing ratio
        Real qv_dry = gas_vmr(i,k,0)*PC::ep_2;
        Real qv_wet = qv_dry/(1.0+qv_dry);
        d_qv(i,k)  = qv_wet;
        // rest of active gases are specified as volume mixing ratios
        d_co2(i,k) = gas_vmr(i,k,1);
        d_o3(i,k)  = gas_vmr(i,k,2);
        d_n2o(i,k) = gas_vmr(i,k,3);
        d_co(i,k)  = gas_vmr(i,k,4);
        d_ch4(i,k) = gas_vmr(i,k,5);
        d_o2(i,k)  = gas_vmr(i,k,6);
        d_n2(i,k)  = gas_vmr(i,k,7);
      });

      d_pint(i,nlay) = p_lev(i,nlay);

      // Compute surface flux from surface temperature
      auto ibot = (p_lay(0,0) > p_lay(0,nlay-1)) ? 0 : nlay;
      d_surf_lw_flux_up(i) = PC::stebol * pow(t_lev(i,ibot), 4.0);
    });
  }
  Kokkos::fence();

  // Run driver
  ad.run(300);

  // Check values; The correct values have been stored in the field manager, we need to
  // copy back to YAKL::Array.
  auto d_sw_flux_up = field_mgr.get_field("sw_flux_up").get_view<Real**>();
  auto d_sw_flux_dn = field_mgr.get_field("sw_flux_dn").get_view<Real**>();
  auto d_sw_flux_dn_dir = field_mgr.get_field("sw_flux_dn_dir").get_view<Real**>();
  auto d_lw_flux_up = field_mgr.get_field("lw_flux_up").get_view<Real**>();
  auto d_lw_flux_dn = field_mgr.get_field("lw_flux_dn").get_view<Real**>();
  auto sw_flux_up_test = real2dk("sw_flux_up_test", ncol, nlay+1);
  auto sw_flux_dn_test = real2dk("sw_flux_dn_test", ncol, nlay+1);
  auto sw_flux_dn_dir_test = real2dk("sw_flux_dn_dir_test",  ncol, nlay+1);
  auto lw_flux_up_test = real2dk("lw_flux_up_test", ncol, nlay+1);
  auto lw_flux_dn_test = real2dk("lw_flux_dn_test", ncol, nlay+1);
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ncol, nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlay+1), [&] (const int& k) {
        if (k < nlay) t_lay(i,k) = d_tmid(i,k);

        sw_flux_up_test(i,k)     = d_sw_flux_up(i,k);
        sw_flux_dn_test(i,k)     = d_sw_flux_dn(i,k);
        sw_flux_dn_dir_test(i,k) = d_sw_flux_dn_dir(i,k);
        lw_flux_up_test(i,k)     = d_lw_flux_up(i,k);
        lw_flux_dn_test(i,k)     = d_lw_flux_dn(i,k);
      });
    });
  }
  Kokkos::fence();

  // Dumb check to verify that we did indeed update temperature
  REQUIRE(Kokkos::create_mirror_view_and_copy(HostDevice(), t_lay)(0,0) !=
          Kokkos::create_mirror_view_and_copy(HostDevice(), T_mid0)(0,0));

  // Make sure fluxes from field manager that were calculated in AD call of RRTMGP match reference fluxes from input file
  // We use all_close here instead of all_equals because we are only able
  // to approximate the solar zenith angle used in the RRTMGP clear-sky
  // test problem with our trial and error lat/lon values, so fluxes will
  // be slightly off. We just verify that they are all "close" here, within
  // some tolerance. Computation of level temperatures from midpoints is
  // also unable to exactly reproduce the values in the test problem, so
  // computed fluxes will be further off from the reference calculation.
  //
  // We need to create local copies with only the columns specific to our rank in case we have split columns over multiple ranks
  auto sw_flux_up_loc     = real2dk("sw_flux_up_loc"    , ncol, nlay+1);
  auto sw_flux_dn_loc     = real2dk("sw_flux_dn_loc"    , ncol, nlay+1);
  auto sw_flux_dn_dir_loc = real2dk("sw_flux_dn_dir_loc", ncol, nlay+1);
  auto lw_flux_up_loc     = real2dk("lw_flux_up_loc"    , ncol, nlay+1);
  auto lw_flux_dn_loc     = real2dk("lw_flux_dn_loc"    , ncol, nlay+1);
  Kokkos::parallel_for(MDRP::template get<2>({nlay+1,ncol}), KOKKOS_LAMBDA(int ilay, int icol) {
    auto icol_all = ncol * irank + icol;
    sw_flux_up_loc(icol,ilay) = sw_flux_up_ref(icol_all,ilay);
    sw_flux_dn_loc(icol,ilay) = sw_flux_dn_ref(icol_all,ilay);
    sw_flux_dn_dir_loc(icol,ilay) = sw_flux_dn_dir_ref(icol_all,ilay);
    lw_flux_up_loc(icol,ilay) = lw_flux_up_ref(icol_all,ilay);
    lw_flux_dn_loc(icol,ilay) = lw_flux_dn_ref(icol_all,ilay);
  });
  REQUIRE(utils_t::all_close(sw_flux_up_loc    , sw_flux_up_test    , 1.0));
  REQUIRE(utils_t::all_close(sw_flux_dn_loc    , sw_flux_dn_test    , 1.0));
  REQUIRE(utils_t::all_close(sw_flux_dn_dir_loc, sw_flux_dn_dir_test, 1.0));
  REQUIRE(utils_t::all_close(lw_flux_up_loc    , lw_flux_up_test    , 1.0));
  REQUIRE(utils_t::all_close(lw_flux_dn_loc    , lw_flux_dn_test    , 1.0));

  // Finalize the driver. YAKL will be finalized inside
  // RRTMGPRadiation::finalize_impl after RRTMGP has had the
  // opportunity to deallocate all it's arrays.
  ad.finalize();
}
#endif

}
