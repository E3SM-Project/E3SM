#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "control/atmosphere_surface_coupling_importer.hpp"
#include "control/atmosphere_surface_coupling_exporter.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/field/field_manager.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include <ekat/ekat_parse_yaml_file.hpp>
#include <ekat/util/ekat_test_utils.hpp>

#include <iomanip>

namespace scream {

using vos_type = std::vector<std::string>;
using vor_type = std::vector<Real>;
constexpr Real test_tol = std::numeric_limits<Real>::epsilon()*1e4;
constexpr Real FillValue = -99999.0;

// Test function for prescribed values
Real test_func(const int col, const int t) {
  return (col + 1 + t);
}

// Wrapper for output manager that also extracts the list of files
class OutputManager4Test : public scream::OutputManager
{
public:
  OutputManager4Test() = default;

  void runme(const util::TimeStamp& ts) {
    run(ts);
    update_file_list();
  }

  std::vector<std::string> get_list_of_files() { return m_list_of_files; }
private:
  void update_file_list() {
    if (std::find(m_list_of_files.begin(),m_list_of_files.end(), m_output_file_specs.filename) == m_list_of_files.end()) {
      m_list_of_files.push_back(m_output_file_specs.filename);
    }
  }
  std::vector<std::string> m_list_of_files;
};

std::vector<std::string> create_from_file_test_data(const ekat::Comm& comm, const util::TimeStamp& t0, const int ncols )
{
  // Create a grids manager on the fly
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_type{"Point Grid"});
  auto& pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_type{"Physics"});
  pl.set("number_of_global_columns",ncols);
  pl.set("number_of_vertical_levels",1); // We don't care about levels for a surface only file
  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();
  // Create a fields manager on the fly with the appropriate fields and grid.
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  const auto grid = gm->get_grid("Physics");
  const int nlcols = grid->get_num_local_dofs();
  const auto dofs_gids = grid->get_dofs_gids().get_view<const int*,Host>();
  std::vector<std::string> fnames = {"lwdn"};
  FieldLayout layout({COL},{nlcols});
  auto fm = std::make_shared<FieldManager>(grid);
  fm->registration_begins();
  fm->registration_ends();
  auto nondim = Units::nondimensional();
  for (auto name : fnames) {
    FieldIdentifier fid(name,layout,nondim,grid->name());
    Field f(fid);
    f.allocate_view();
    // Initialize data
    auto f_view_h = f.get_view<Real*,Host>();
    for (int ii=0; ii<nlcols; ii++) {
      f_view_h(ii) = test_func(ii,0);
    }
    f.sync_to_dev();
    // Update timestamp
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  // Create output manager to handle the data
  scorpio::init_subsystem(comm);
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix",std::string("surface_coupling_forcing"));
  om_pl.set("Field Names",fnames);
  om_pl.set("Averaging Type", std::string("INSTANT"));
  om_pl.set("Max Snapshots Per File",2);
  om_pl.set<double>("fill_value",FillValue);
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units",std::string("nsteps"));
  ctrl_pl.set("Frequency",1);
  ctrl_pl.set("save_grid_data",false);
  OutputManager4Test om;
  om.initialize(comm,om_pl,t0,false);
  om.setup(fm,gm);
  // Create output data:
  // T=3600, well above the max timestep for the test.
  auto tw = t0;
  const int dt = 3600;
  for (auto name : fnames) {
    auto field = fm->get_field(name);
    // Note we only care about surface values so we only need to generate data over nlcols.
    auto f_view_h = field.get_view<Real*,Host>();
    for (int ii=0; ii<nlcols; ii++) {
      f_view_h(ii) = test_func(ii,dt);
    }
    field.sync_to_dev();
  }
  tw+=dt;
  om.runme(tw);

  // Done, finalize
  om.finalize();
  scorpio::finalize_subsystem();
  return om.get_list_of_files();
}

void setup_import_and_export_data(
        std::mt19937_64& engine,
  const ekat::Comm& comm,
  // Imports
  const int num_cpl_imports, const int num_scream_imports,
  const KokkosTypes<HostDevice>::view_1d<int >& import_cpl_indices_view,
  const KokkosTypes<HostDevice>::view_1d<int >& import_vec_comps_view,
  const KokkosTypes<HostDevice>::view_1d<Real>& import_constant_multiple_view,
  const KokkosTypes<HostDevice>::view_1d<bool>& do_import_during_init_view,
  // Export data
  const int num_cpl_exports, const int num_scream_exports,
  const KokkosTypes<HostDevice>::view_1d<int >& export_cpl_indices_view,
  const KokkosTypes<HostDevice>::view_1d<int >& export_vec_comps_view,
  const KokkosTypes<HostDevice>::view_1d<Real>& export_constant_multiple_view,
  const KokkosTypes<HostDevice>::view_1d<bool>& do_export_during_init_view)
{
  // IMPORTS

  // Default values
  Kokkos::deep_copy(import_vec_comps_view,         -1);
  Kokkos::deep_copy(import_constant_multiple_view,  1);
  Kokkos::deep_copy(do_import_during_init_view,     false);

  // Set cpl indices. We do a random shuffle of [0, num_cpl_imports),
  // then assign the scream imports to the first num_scream_imports
  // in the shuffled indices.
  {
    std::vector<int> import_order(num_cpl_imports);
    for (int f=0; f<num_cpl_imports; ++f) { import_order[f] = f; }
    std::shuffle(import_order.begin(), import_order.end(),engine);

    comm.broadcast(import_order.data(),num_cpl_imports,0);
    for (int f=0; f<num_scream_imports; ++f) {
      import_cpl_indices_view(f) = import_order[f];
    }
  }

  // Set vector components
  import_vec_comps_view(10) = 0;
  import_vec_comps_view(11) = 1;

  // Set constant multiples
  import_constant_multiple_view(9)  = -1;
  import_constant_multiple_view(10) = -1;
  import_constant_multiple_view(11) = -1;
  import_constant_multiple_view(12) = -1;
  import_constant_multiple_view(13) = -1;

  // Set boolean for importing during intialization
  do_import_during_init_view(10) = true;
  do_import_during_init_view(11) = true;
  do_import_during_init_view(12) = true;
  do_import_during_init_view(13) = true;

  // EXPORTS

  // Default values
  Kokkos::deep_copy(export_vec_comps_view,         -1);
  Kokkos::deep_copy(export_constant_multiple_view,  1);
  Kokkos::deep_copy(do_export_during_init_view,     false);

  // Set cpl indices. We do a random shuffle of [0, num_cpl_exports),
  // then assign the scream exports to the first num_scream_exports
  // in the shuffled indices.
  {
    std::vector<int> export_order(num_cpl_exports);
    for (int f=0; f<num_cpl_exports; ++f) { export_order[f] = f; }
    std::shuffle(export_order.begin(), export_order.end(),engine);
    comm.broadcast(export_order.data(),num_cpl_exports,0);
    for (int f=0; f<num_scream_exports; ++f) {
      export_cpl_indices_view(f) = export_order[f];
    }
  }

  // Set boolean for exporting during intialization
  do_export_during_init_view(0) = true;
  do_export_during_init_view(1) = true;
  do_export_during_init_view(2) = true;
  do_export_during_init_view(3) = true;
  do_export_during_init_view(4) = true;
  do_export_during_init_view(5) = true;
  do_export_during_init_view(6) = true;
  do_export_during_init_view(7) = true;
  do_export_during_init_view(8) = true;
}

void test_imports(const FieldManager& fm,
                  const KokkosTypes<HostDevice>::view_2d<Real> import_data_view,
                  const KokkosTypes<HostDevice>::view_1d<int>  import_cpl_indices_view,
                  const KokkosTypes<HostDevice>::view_1d<Real> import_constant_multiple_view,
                  const bool called_directly_after_init = false)
{
  // Sync to host for comparing to import data
  fm.get_field("sfc_alb_dir_vis" ).sync_to_host();
  fm.get_field("sfc_alb_dir_nir" ).sync_to_host();
  fm.get_field("sfc_alb_dif_vis" ).sync_to_host();
  fm.get_field("sfc_alb_dif_nir" ).sync_to_host();
  fm.get_field("surf_radiative_T").sync_to_host();
  fm.get_field("T_2m"            ).sync_to_host();
  fm.get_field("qv_2m"           ).sync_to_host();
  fm.get_field("wind_speed_10m"  ).sync_to_host();
  fm.get_field("snow_depth_land" ).sync_to_host();
  fm.get_field("surf_lw_flux_up" ).sync_to_host();
  fm.get_field("surf_mom_flux"   ).sync_to_host();
  fm.get_field("surf_sens_flux"  ).sync_to_host();
  fm.get_field("surf_evap"       ).sync_to_host();
  fm.get_field("ocnfrac"         ).sync_to_host();
  fm.get_field("landfrac"        ).sync_to_host();
  const auto sfc_alb_dir_vis  = fm.get_field("sfc_alb_dir_vis" ).get_view<const Real*,  Host>();
  const auto sfc_alb_dir_nir  = fm.get_field("sfc_alb_dir_nir" ).get_view<const Real*,  Host>();
  const auto sfc_alb_dif_vis  = fm.get_field("sfc_alb_dif_vis" ).get_view<const Real*,  Host>();
  const auto sfc_alb_dif_nir  = fm.get_field("sfc_alb_dif_nir" ).get_view<const Real*,  Host>();
  const auto surf_radiative_T = fm.get_field("surf_radiative_T").get_view<const Real*,  Host>();
  const auto T_2m             = fm.get_field("T_2m"            ).get_view<const Real*,  Host>();
  const auto qv_2m            = fm.get_field("qv_2m"           ).get_view<const Real*,  Host>();
  const auto wind_speed_10m   = fm.get_field("wind_speed_10m"  ).get_view<const Real*,  Host>();
  const auto snow_depth_land  = fm.get_field("snow_depth_land" ).get_view<const Real*,  Host>();
  const auto surf_lw_flux_up  = fm.get_field("surf_lw_flux_up" ).get_view<const Real*,  Host>();
  const auto surf_mom_flux    = fm.get_field("surf_mom_flux"   ).get_view<const Real**, Host>();
  const auto surf_sens_flux   = fm.get_field("surf_sens_flux"  ).get_view<const Real*,  Host>();
  const auto surf_evap        = fm.get_field("surf_evap"       ).get_view<const Real*,  Host>();
  const auto ocnfrac          = fm.get_field("ocnfrac"         ).get_view<const Real*,  Host>();
  const auto landfrac         = fm.get_field("landfrac"        ).get_view<const Real*,  Host>();

  const int ncols = surf_evap.extent(0);

  for (int i=0; i<ncols; ++i) {

    // Check cpl data to scream fields

    // The following are imported both during initialization and run phase
    EKAT_REQUIRE(surf_mom_flux(i,0) == import_constant_multiple_view(10)*import_data_view(i, import_cpl_indices_view(10)));
    EKAT_REQUIRE(surf_mom_flux(i,1) == import_constant_multiple_view(11)*import_data_view(i, import_cpl_indices_view(11)));
    EKAT_REQUIRE(surf_sens_flux(i)  == import_constant_multiple_view(12)*import_data_view(i, import_cpl_indices_view(12)));
    EKAT_REQUIRE(surf_evap(i)       == import_constant_multiple_view(13)*import_data_view(i, import_cpl_indices_view(13)));

    // The following are only imported during run phase. If this test is called
    // during initialization, all values should be the default value.
    // TODO: Why are some of these FillValue and others are 0.0?  For the former,
    // they are gathered from output it seems, so they take the FillValue.  While
    // the few ones with 0.0 seem to take there initial value from the field initialization
    // which is 0.0 I believe.  Still, based on the comment none of these should be read in
    // yet right?
    if (called_directly_after_init) {
      EKAT_REQUIRE(sfc_alb_dir_vis(i)  == FillValue);
      EKAT_REQUIRE(sfc_alb_dir_nir(i)  == FillValue);
      EKAT_REQUIRE(sfc_alb_dif_vis(i)  == FillValue);
      EKAT_REQUIRE(sfc_alb_dif_nir(i)  == FillValue);
      EKAT_REQUIRE(surf_radiative_T(i) == 0.0);
      EKAT_REQUIRE(T_2m(i)             == 0.0);
      EKAT_REQUIRE(qv_2m(i)            == 0.0);
      EKAT_REQUIRE(wind_speed_10m(i)   == 0.0);
      EKAT_REQUIRE(snow_depth_land(i)  == 0.0);
      EKAT_REQUIRE(surf_lw_flux_up(i)  == FillValue);
      EKAT_REQUIRE(ocnfrac(i)          == FillValue);
      EKAT_REQUIRE(landfrac(i)         == FillValue);
    } else {
      EKAT_REQUIRE(sfc_alb_dir_vis(i)  == import_constant_multiple_view(0 )*import_data_view(i, import_cpl_indices_view(0)));
      EKAT_REQUIRE(sfc_alb_dir_nir(i)  == import_constant_multiple_view(1 )*import_data_view(i, import_cpl_indices_view(1)));
      EKAT_REQUIRE(sfc_alb_dif_vis(i)  == import_constant_multiple_view(2 )*import_data_view(i, import_cpl_indices_view(2)));
      EKAT_REQUIRE(sfc_alb_dif_nir(i)  == import_constant_multiple_view(3 )*import_data_view(i, import_cpl_indices_view(3)));
      EKAT_REQUIRE(surf_radiative_T(i) == import_constant_multiple_view(4 )*import_data_view(i, import_cpl_indices_view(4)));
      EKAT_REQUIRE(T_2m(i)             == import_constant_multiple_view(5 )*import_data_view(i, import_cpl_indices_view(5)));
      EKAT_REQUIRE(qv_2m(i)            == import_constant_multiple_view(6 )*import_data_view(i, import_cpl_indices_view(6)));
      EKAT_REQUIRE(wind_speed_10m(i)   == import_constant_multiple_view(7 )*import_data_view(i, import_cpl_indices_view(7)));
      EKAT_REQUIRE(snow_depth_land(i)  == import_constant_multiple_view(8 )*import_data_view(i, import_cpl_indices_view(8)));
      EKAT_REQUIRE(surf_lw_flux_up(i)  == import_constant_multiple_view(9 )*import_data_view(i, import_cpl_indices_view(9)));
      EKAT_REQUIRE(ocnfrac(i)          == import_constant_multiple_view(14)*import_data_view(i, import_cpl_indices_view(14)));
      EKAT_REQUIRE(landfrac(i)         == import_constant_multiple_view(15)*import_data_view(i, import_cpl_indices_view(15)));
    }
  }
}

void test_exports(const FieldManager& fm,
                  const KokkosTypes<HostDevice>::view_2d<Real> export_data_view,
                  const KokkosTypes<HostDevice>::view_1d<int>  export_cpl_indices_view,
                  const KokkosTypes<HostDevice>::view_1d<Real> export_constant_multiple_view,
                  const ekat::ParameterList prescribed_constants,
                  const int dt,
                  const bool called_directly_after_init = false)
{
  using PF = PhysicsFunctions<DefaultDevice>;
  using PC = physics::Constants<Real>;

  // Some computed fields rely on calculations that are done in the AD.
  // Recompute here and verify that they were exported correctly.
  const auto pseudo_density       = fm.get_field("pseudo_density").get_view<const Real**>();
  const auto p_mid                = fm.get_field("p_mid").get_view<const Real**>();
  const auto p_int                = fm.get_field("p_int").get_view<const Real**>();
  const auto T_mid                = fm.get_field("T_mid").get_view<const Real**>();
  const auto qv                   = fm.get_field("qv").get_view<const Real**>();
  const auto phis                 = fm.get_field("phis").get_view<const Real*>();
  const auto precip_liq_surf_mass = fm.get_field("precip_liq_surf_mass").get_view<const Real*>();
  const auto precip_ice_surf_mass = fm.get_field("precip_ice_surf_mass").get_view<const Real*>();

  const int ncols = fm.get_grid()->get_num_local_dofs();
  const int nlevs = fm.get_grid()->get_num_vertical_levels();

  KokkosTypes<DefaultDevice>::view_2d<Real> dz   ("dz   ", ncols, nlevs);
  KokkosTypes<DefaultDevice>::view_2d<Real> z_mid("z_mid", ncols, nlevs);
  KokkosTypes<DefaultDevice>::view_2d<Real> z_int("z_int", ncols, nlevs+1);

  KokkosTypes<DefaultDevice>::view_1d<Real> Sa_z      ("Sa_z",       ncols);
  KokkosTypes<DefaultDevice>::view_1d<Real> Sa_ptem   ("Sa_ptem",    ncols);
  KokkosTypes<DefaultDevice>::view_1d<Real> Sa_dens   ("Sa_dens",    ncols);
  KokkosTypes<DefaultDevice>::view_1d<Real> Sa_pslv   ("Sa_pslv",    ncols);
  KokkosTypes<DefaultDevice>::view_1d<Real> Faxa_rainl("Faxa_rainl", ncols);
  KokkosTypes<DefaultDevice>::view_1d<Real> Faxa_snowl("Faxa_snowl", ncols);

  const auto setup_policy =
      ekat::ExeSpaceUtils<KokkosTypes<DefaultDevice>::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncols, nlevs);
  Kokkos::parallel_for(setup_policy, KOKKOS_LAMBDA(const Kokkos::TeamPolicy<KokkosTypes<DefaultDevice>::ExeSpace>::member_type& team) {
    const int i = team.league_rank();

    const auto pseudo_density_i  = ekat::subview(pseudo_density, i);
    const auto p_mid_i           = ekat::subview(p_mid,          i);
    const auto p_int_i           = ekat::subview(p_int,          i);
    const auto T_mid_i           = ekat::subview(T_mid,          i);
    const auto qv_i              = ekat::subview(qv,             i);
    const auto dz_i              = ekat::subview(dz,             i);
    const auto z_int_i           = ekat::subview(z_int,          i);
    const auto z_mid_i           = ekat::subview(z_mid,          i);

    // Compute vertical layer thickness
    PF::calculate_dz(team, pseudo_density_i, p_mid_i, T_mid_i, qv_i, dz_i);
    team.team_barrier();

    // Compute vertical layer heights (relative to ground surface rather than from sea level).
    // Use z_int(nlevs) = z_surf = 0.0.
    const Real z_surf = 0.0;
    PF::calculate_z_int(team, nlevs, dz_i, z_surf, z_int_i);
    team.team_barrier();
    PF::calculate_z_mid(team, nlevs, z_int_i, z_mid_i);
    team.team_barrier();

    // Calculate air temperature at bottom of cell closest to the ground for PSL
    const Real T_int_bot = PF::calculate_surface_air_T(T_mid_i(nlevs-1),z_mid_i(nlevs-1));

    Sa_z(i)       = z_mid_i(nlevs-1);
    Sa_ptem(i)    = PF::calculate_theta_from_T(T_mid_i(nlevs-1), p_mid_i(nlevs-1));
    Sa_dens(i)    = PF::calculate_density(pseudo_density_i(nlevs-1), dz_i(nlevs-1));
    Sa_pslv(i)    = PF::calculate_psl(T_int_bot, p_int_i(nlevs), phis(i));

    if (not called_directly_after_init) {
      Faxa_rainl(i) = precip_liq_surf_mass(i)/dt*(1000.0/PC::RHO_H2O);
      Faxa_snowl(i) = precip_ice_surf_mass(i)/dt*(1000.0/PC::RHO_H2O);
    }
  });

  // Sync to host for comparing to export data
  fm.get_field("p_mid"           ).sync_to_host();
  fm.get_field("T_mid"           ).sync_to_host();
  fm.get_field("qv"              ).sync_to_host();
  fm.get_field("horiz_winds"     ).sync_to_host();
  fm.get_field("sfc_flux_dir_nir").sync_to_host();
  fm.get_field("sfc_flux_dir_vis").sync_to_host();
  fm.get_field("sfc_flux_dif_nir").sync_to_host();
  fm.get_field("sfc_flux_dif_vis").sync_to_host();
  fm.get_field("sfc_flux_sw_net" ).sync_to_host();
  fm.get_field("sfc_flux_lw_dn"  ).sync_to_host();
  const auto p_mid_h            = fm.get_field("p_mid"           ).get_view<const Real**,  Host>();
  const auto T_mid_h            = fm.get_field("T_mid"           ).get_view<const Real**,  Host>();
  const auto qv_h               = fm.get_field("qv"              ).get_view<const Real**,  Host>();
  const auto horiz_winds_h      = fm.get_field("horiz_winds"     ).get_view<const Real***, Host>();
  const auto sfc_flux_dir_nir_h = fm.get_field("sfc_flux_dir_nir").get_view<const Real*,   Host>();
  const auto sfc_flux_dir_vis_h = fm.get_field("sfc_flux_dir_vis").get_view<const Real*,   Host>();
  const auto sfc_flux_dif_nir_h = fm.get_field("sfc_flux_dif_nir").get_view<const Real*,   Host>();
  const auto sfc_flux_dif_vis_h = fm.get_field("sfc_flux_dif_vis").get_view<const Real*,   Host>();
  const auto sfc_flux_sw_net_h  = fm.get_field("sfc_flux_sw_net" ).get_view<const Real*,   Host>();
  const auto sfc_flux_lw_dn_h   = fm.get_field("sfc_flux_lw_dn"  ).get_view<const Real*,   Host>();
  const auto Sa_z_h       = Kokkos::create_mirror_view_and_copy(HostDevice(), Sa_z);
  const auto Sa_ptem_h    = Kokkos::create_mirror_view_and_copy(HostDevice(), Sa_ptem);
  const auto Sa_dens_h    = Kokkos::create_mirror_view_and_copy(HostDevice(), Sa_dens);
  const auto Sa_pslv_h    = Kokkos::create_mirror_view_and_copy(HostDevice(), Sa_pslv);
  const auto Faxa_rainl_h = Kokkos::create_mirror_view_and_copy(HostDevice(), Faxa_rainl);
  const auto Faxa_snowl_h = Kokkos::create_mirror_view_and_copy(HostDevice(), Faxa_snowl);

  // Recall that two fields have been set to export to a constant value, so we load those constants from the parameter list here:
  using vor_type = std::vector<Real>;
  const auto prescribed_const_values = prescribed_constants.get<vor_type>("values");
  const Real Faxa_swndf_const = prescribed_const_values[0];
  const Real Faxa_swndv_const = prescribed_const_values[1];


  // Check cpl data to scream fields
  for (int i=0; i<ncols; ++i) {
    const Real Faxa_lwdn_file = test_func(i,dt);

    // The following are exported both during initialization and run phase
    EKAT_REQUIRE(export_constant_multiple_view(0)*Sa_z_h(i)                  == export_data_view(i, export_cpl_indices_view(0)));
    EKAT_REQUIRE(export_constant_multiple_view(1)*horiz_winds_h(i,0,nlevs-1) == export_data_view(i, export_cpl_indices_view(1)));
    EKAT_REQUIRE(export_constant_multiple_view(2)*horiz_winds_h(i,1,nlevs-1) == export_data_view(i, export_cpl_indices_view(2)));
    EKAT_REQUIRE(export_constant_multiple_view(3)*T_mid_h(i,nlevs-1)         == export_data_view(i, export_cpl_indices_view(3)));
    EKAT_REQUIRE(export_constant_multiple_view(4)*Sa_ptem_h(i)               == export_data_view(i, export_cpl_indices_view(4)));
    EKAT_REQUIRE(export_constant_multiple_view(5)*p_mid_h(i,nlevs-1)         == export_data_view(i, export_cpl_indices_view(5)));
    EKAT_REQUIRE(export_constant_multiple_view(6)*qv_h(i,nlevs-1)            == export_data_view(i, export_cpl_indices_view(6)));
    EKAT_REQUIRE(export_constant_multiple_view(7)*Sa_dens_h(i)               == export_data_view(i, export_cpl_indices_view(7)));
    EKAT_REQUIRE(export_constant_multiple_view(8)*Sa_pslv_h(i)               == export_data_view(i, export_cpl_indices_view(8)));

    // The following are only exported during run phase. If this test is called
    // during initialization, all values should have been set to 0.
    if (called_directly_after_init) {
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(9 )));
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(10)));
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(11)));
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(12)));
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(13)));
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(14)));
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(15)));
      EKAT_REQUIRE(0 == export_data_view(i, export_cpl_indices_view(16)));
    } else {
      EKAT_REQUIRE(export_constant_multiple_view(9 )*Faxa_rainl_h(i)       == export_data_view(i, export_cpl_indices_view(9 )));
      EKAT_REQUIRE(export_constant_multiple_view(10)*Faxa_snowl_h(i)       == export_data_view(i, export_cpl_indices_view(10)));
      EKAT_REQUIRE(export_constant_multiple_view(11)*sfc_flux_dir_nir_h(i) == export_data_view(i, export_cpl_indices_view(11)));
      EKAT_REQUIRE(export_constant_multiple_view(12)*sfc_flux_dir_vis_h(i) == export_data_view(i, export_cpl_indices_view(12)));
      EKAT_REQUIRE(Faxa_swndf_const                                        == export_data_view(i, export_cpl_indices_view(13)));
      EKAT_REQUIRE(Faxa_swndv_const                                        == export_data_view(i, export_cpl_indices_view(14)));
      EKAT_REQUIRE(export_constant_multiple_view(15)*sfc_flux_sw_net_h(i)  == export_data_view(i, export_cpl_indices_view(15)));
      EKAT_REQUIRE(std::abs(Faxa_lwdn_file - export_data_view(i, export_cpl_indices_view(16)))<test_tol);
    }
  }
}

TEST_CASE("surface-coupling", "") {
  using namespace scream;
  using namespace scream::control;

  // Create a comm
  ekat::Comm atm_comm (MPI_COMM_WORLD);
  auto engine = setup_random_test(&atm_comm);

  // Load ad parameter list
  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname,ad_params);

  // Parameters
  auto& ts          = ad_params.sublist("time_stepping");
  const auto t0_str = ts.get<std::string>("run_t0");
  const auto t0     = util::str_to_time_stamp(t0_str);
  const auto gmp    = ad_params.sublist("grids_manager");
  const auto grid_name = gmp.get<vos_type>("grids_names");
  const auto gdp    = gmp.sublist(grid_name[0]);
  const auto ncol_in = gdp.get<int>("number_of_global_columns");

  // Set two export fields to be randomly set to a constant
  // This requires us to add a sublist to the parsed AD params yaml list.
  std::uniform_real_distribution<Real> pdf_real_constant_data(0.0,1.0);

  auto& ap_params     = ad_params.sublist("atmosphere_processes");
  auto& sc_exp_params = ap_params.sublist("SurfaceCouplingExporter");
  // Set up forcing to a constant value
  const Real Faxa_swndf_const = pdf_real_constant_data(engine);
  const Real Faxa_swvdf_const = pdf_real_constant_data(engine);
  const vos_type exp_const_fields = {"Faxa_swndf","Faxa_swvdf"};
  vor_type exp_const_values = {Faxa_swndf_const,Faxa_swvdf_const};
  atm_comm.broadcast(exp_const_values.data(),2,0);
  auto& exp_const_params = sc_exp_params.sublist("prescribed_constants");
  exp_const_params.set<vos_type>("fields",exp_const_fields);
  exp_const_params.set<vor_type>("values",exp_const_values);
  // Set up forcing to data interpolated from file
  const auto exp_file_files = create_from_file_test_data(atm_comm, t0, ncol_in);
  const vos_type exp_file_fields = {"Faxa_lwdn"};
  // Test the use of an alternative name as stored in the data file(s).
  const vos_type exp_file_fields_alt_name = {"Faxa_lwdn:lwdn"};
  auto& exp_file_params = sc_exp_params.sublist("prescribed_from_file");
  exp_file_params.set<vos_type>("fields",exp_file_fields);
  exp_file_params.set<vos_type>("fields_alt_name",exp_file_fields_alt_name);
  exp_file_params.set<vos_type>("files",exp_file_files);

  // Need to register products in the factory *before* we create any atm process or grids manager.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("SurfaceCouplingImporter",&create_atmosphere_process<SurfaceCouplingImporter>);
  proc_factory.register_product("SurfaceCouplingExporter",&create_atmosphere_process<SurfaceCouplingExporter>);
  register_mesh_free_grids_manager();
  register_diagnostics();

  // Create the AD
  AtmosphereDriver ad;
  ad.set_comm(atm_comm);
  ad.set_params(ad_params);
  ad.init_scorpio ();
  ad.init_time_stamps (t0, t0);
  ad.create_output_managers ();
  ad.create_atm_processes ();
  ad.create_grids ();
  ad.create_fields ();

  const int   ncols = ad.get_grids_manager()->get_grid("Physics")->get_num_local_dofs();

  // Create engine and pdfs for random test data
  std::uniform_int_distribution<int> pdf_int_additional_fields(0,10);
  std::uniform_int_distribution<int> pdf_int_dt(1,1800);
  std::uniform_real_distribution<Real> pdf_real_import_data(0.0,1.0);
  // Set up random value for dt
  int dt = pdf_int_dt(engine);
  atm_comm.broadcast(&dt,1,0);
  // Setup views to test import/export. For this test we consider a random number of non-imported/exported
  // cpl fields (in addition to the required scream imports/exports), then assign a random, non-repeating
  // cpl index for each field in [0, num_cpl_fields).
  const int num_scream_imports = 16;
  const int num_scream_exports = 17;
  KokkosTypes<HostDevice>::view_1d<int> additional_import_exports("additional_import_exports", 2);
  ekat::genRandArray(additional_import_exports, engine, pdf_int_additional_fields);
  atm_comm.broadcast(additional_import_exports.data(),2,0);
  const int num_additional_imports = additional_import_exports(0);
  const int num_additional_exports = additional_import_exports(1);
  const int num_cpl_imports = num_scream_imports + num_additional_imports;
  const int num_cpl_exports = num_scream_exports + num_additional_exports;
  // Import data is of size num_cpl_imports, the rest of the views are size num_scream_imports.
  KokkosTypes<HostDevice>::view_2d<Real> import_data_view             ("import_data",
                                                                       ncols, num_cpl_imports);
  KokkosTypes<HostDevice>::view_1d<int>  import_cpl_indices_view      ("import_cpl_indices",
                                                                       num_scream_imports);
  KokkosTypes<HostDevice>::view_1d<int>  import_vec_comps_view        ("import_vec_comps",
                                                                       num_scream_imports);
  KokkosTypes<HostDevice>::view_1d<Real> import_constant_multiple_view("import_constant_multiple_view",
                                                                       num_scream_imports);
  KokkosTypes<HostDevice>::view_1d<bool> do_import_during_init_view   ("do_import_during_init_view",
                                                                       num_scream_imports);
  // Set import data to random (0,1) values
  ekat::genRandArray(import_data_view, engine, pdf_real_import_data);
  atm_comm.broadcast(import_data_view.data(), num_cpl_imports, 0);
  // Set import names
  char import_names[num_scream_imports][32];
  std::strcpy(import_names[0],  "sfc_alb_dir_vis");
  std::strcpy(import_names[1],  "sfc_alb_dir_nir");
  std::strcpy(import_names[2],  "sfc_alb_dif_vis");
  std::strcpy(import_names[3],  "sfc_alb_dif_nir");
  std::strcpy(import_names[4],  "surf_radiative_T");
  std::strcpy(import_names[5],  "T_2m");
  std::strcpy(import_names[6],  "qv_2m");
  std::strcpy(import_names[7],  "wind_speed_10m");
  std::strcpy(import_names[8],  "snow_depth_land");
  std::strcpy(import_names[9],  "surf_lw_flux_up");
  std::strcpy(import_names[10], "surf_mom_flux");
  std::strcpy(import_names[11], "surf_mom_flux");
  std::strcpy(import_names[12], "surf_sens_flux");
  std::strcpy(import_names[13], "surf_evap");
  std::strcpy(import_names[14], "ocnfrac");
  std::strcpy(import_names[15], "landfrac");

  // Export data is of size num_cpl_exports, the rest of the views are size num_scream_exports.
  KokkosTypes<HostDevice>::view_2d<Real> export_data_view             ("export_data",
                                                                       ncols, num_cpl_exports);
  KokkosTypes<HostDevice>::view_1d<int>  export_cpl_indices_view      ("export_vec_comps",
                                                                       num_scream_exports);
  KokkosTypes<HostDevice>::view_1d<int>  export_vec_comps_view        ("export_vec_comps",
                                                                       num_scream_exports);
  KokkosTypes<HostDevice>::view_1d<Real> export_constant_multiple_view("export_constant_multiple_view",
                                                                       num_scream_exports);
  KokkosTypes<HostDevice>::view_1d<bool> do_export_during_init_view   ("do_export_during_init_view",
                                                                       num_scream_exports);
  // Set export data to -1. All values should be overwritten after the run phase.
  Kokkos::deep_copy(export_data_view, -1.0);
  // Set names. For all non-scream exports, set to 0.
  char export_names[num_scream_exports][32];
  std::strcpy(export_names[0],  "Sa_z"       );
  std::strcpy(export_names[1],  "Sa_u"       );
  std::strcpy(export_names[2],  "Sa_v"       );
  std::strcpy(export_names[3],  "Sa_tbot"    );
  std::strcpy(export_names[4],  "Sa_ptem"    );
  std::strcpy(export_names[5],  "Sa_pbot"    );
  std::strcpy(export_names[6],  "Sa_shum"    );
  std::strcpy(export_names[7],  "Sa_dens"    );
  std::strcpy(export_names[8],  "Sa_pslv"    );
  std::strcpy(export_names[9],  "Faxa_rainl" );
  std::strcpy(export_names[10], "Faxa_snowl" );
  std::strcpy(export_names[11], "Faxa_swndr" );
  std::strcpy(export_names[12], "Faxa_swvdr" );
  std::strcpy(export_names[13], "Faxa_swndf" );
  std::strcpy(export_names[14], "Faxa_swvdf" );
  std::strcpy(export_names[15], "Faxa_swnet" );
  std::strcpy(export_names[16], "Faxa_lwdn"  );

  // Setup the import/export data. This is meant to replicate the structures coming
  // from mct_coupling/scream_cpl_indices.F90
  setup_import_and_export_data(engine, atm_comm,
                               num_cpl_imports, num_scream_imports,
                               import_cpl_indices_view, import_vec_comps_view,
                               import_constant_multiple_view, do_import_during_init_view,
                               num_cpl_exports, num_scream_exports,
                               export_cpl_indices_view, export_vec_comps_view,
                               export_constant_multiple_view, do_export_during_init_view);

  // Have the AD setup surface coupling data
  ad.setup_surface_coupling_data_manager(SurfaceCouplingTransferType::Import,
                                         num_cpl_imports, num_scream_imports, ncols, import_data_view.data(),
                                         import_names[0], import_cpl_indices_view.data(), import_vec_comps_view.data(),
                                         import_constant_multiple_view.data(), do_import_during_init_view.data());
  ad.setup_surface_coupling_data_manager(SurfaceCouplingTransferType::Export,
                                         num_cpl_exports, num_scream_exports, ncols, export_data_view.data(),
                                         export_names[0], export_cpl_indices_view.data(), export_vec_comps_view.data(),
                                         export_constant_multiple_view.data(), do_export_during_init_view.data());

  // Initialize the AD
  ad.initialize_fields ();
  ad.initialize_output_managers ();
  ad.initialize_atm_procs ();

  const auto fm = ad.get_field_mgr("Physics");

  // Verify any initial imports/exports were done as expected
  test_imports(*fm, import_data_view, import_cpl_indices_view,
               import_constant_multiple_view, true);
  test_exports(*fm, export_data_view, export_cpl_indices_view,
               export_constant_multiple_view,  exp_const_params, dt, true);

  // Run the AD
  ad.run(dt);

  // Verify all imports/exports were done as expected
  test_imports(*fm, import_data_view, import_cpl_indices_view,
               import_constant_multiple_view);
  test_exports(*fm, export_data_view, export_cpl_indices_view,
               export_constant_multiple_view, exp_const_params, dt);

  // Finalize  the AD
  ad.finalize();
}

} // empty namespace
