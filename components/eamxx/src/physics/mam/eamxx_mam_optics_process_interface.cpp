#include <ekat/ekat_assert.hpp>
#include <physics/mam/eamxx_mam_optics_process_interface.hpp>
#include <share/property_checks/field_lower_bound_check.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

#include "eamxx_config.h"  // for SCREAM_CIME_BUILD
#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

namespace scream {

MAMOptics::MAMOptics(const ekat::Comm &comm, const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params), aero_config_() {
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);
}

std::string MAMOptics::name() const { return "mam4_optics"; }

void MAMOptics::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();
  Units n_unit(1 / kg, "#/kg");  // number mixing ratios [# / kg air]
  const auto m2 = pow(m, 2);
  const auto s2 = pow(s, 2);

  ncol_     = grid_->get_num_local_dofs();  // number of columns on this rank
  nlev_     = grid_->get_num_vertical_levels();  // number of levels per column
  nswbands_ = mam4::modal_aer_opt::nswbands;     // number of shortwave bands
  nlwbands_ = mam4::modal_aer_opt::nlwbands;     // number of longwave bands

  // Define the different field layouts that will be used for this process

  // Define aerosol optics fields computed by this process.
  auto nondim = Units::nondimensional();
  // 3D layout for short/longwave aerosol fields: columns, number of
  // short/longwave band, nlev
  FieldLayout scalar3d_swband =
      grid_->get_3d_vector_layout(true, nswbands_, "swband");
  FieldLayout scalar3d_lwband =
      grid_->get_3d_vector_layout(true, nlwbands_, "lwband");

  // layout for 3D (2d horiz X 1d vertical) variables at level
  // midpoints/interfaces
  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = grid_->get_3d_scalar_layout(false);
  add_tracers_wet_atm();
  add_fields_dry_atm();

  // layout for 2D (1d horiz X 1d vertical) variables
  FieldLayout scalar2d = grid_->get_2d_scalar_layout();

  add_field<Required>("pseudo_density_dry", scalar3d_mid, Pa, grid_name);

  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  // shortwave aerosol scattering asymmetry parameter [unitless]
  add_field<Computed>("aero_g_sw", scalar3d_swband, nondim, grid_name);

  // shortwave aerosol single-scattering albedo [unitless]
  add_field<Computed>("aero_ssa_sw", scalar3d_swband, nondim, grid_name);

  // shortwave aerosol extinction optical depth
  add_field<Computed>("aero_tau_sw", scalar3d_swband, nondim, grid_name);

  // longwave aerosol extinction optical depth [unitless]
  add_field<Computed>("aero_tau_lw", scalar3d_lwband, nondim, grid_name);

  add_field<Computed>("aodvis", scalar2d, nondim, grid_name);

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  // add tracers, e.g., num_a1, soa_a1
  add_tracers_interstitial_aerosol();
  // add tracer gases, e.g., O3
  add_tracers_gases();
  // add fields e.g., num_c1, soa_c1
  add_fields_cloudborne_aerosol();

  // aerosol-related gases: mass mixing ratios
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_tracer<Updated>(gas_mmr_field_name, grid_, kg / kg);
  }
}

size_t MAMOptics::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

void MAMOptics::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMMOptics.");
}

void MAMOptics::initialize_impl(const RunType run_type) {
  // Check the interval values for the following fields used by this interface.
  // NOTE: We do not include aerosol and gas species, e.g., soa_a1, num_a1,
  // because we automatically added these fields.
  const std::map<std::string, std::pair<Real, Real>> ranges_optics = {
      // optics
      {"pseudo_density_dry", {-1e10, 1e10}},  // FIXME
      {"aero_g_sw", {-1e10, 1e10}},           // FIXME
      {"aero_ssa_sw", {-1e10, 1e10}},         // FIXME
      {"aero_tau_lw", {-1e10, 1e10}},         // FIXME
      {"aero_tau_sw", {-1e10, 1e10}},         // FIXME
      {"aodvis", {-1e10, 1e10}}               // FIXME
  };
  set_ranges_process(ranges_optics);
  add_interval_checks();
  // populate the wet and dry atmosphere states with views from fields and
  // the buffer
  constexpr int ntot_amode = mam4::AeroConfig::num_modes();

  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);
  p_del_ = get_field_in("pseudo_density").get_view<const Real **>();

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  // It populates wet_aero struct (wet_aero_) with:
  // interstitial aerosol, e.g., soa_a_1
  populate_interstitial_wet_aero(wet_aero_);
  // gases, e.g., O3
  populate_gases_wet_aero(wet_aero_);
  // cloudborne aerosol, e.g., soa_c_1
  populate_cloudborne_wet_aero(wet_aero_);
  // It populates dry_aero struct (dry_aero_) with:
  // interstitial aerosol, e.g., soa_a_1
  populate_interstitial_dry_aero(dry_aero_, buffer_);
  // gases, e.g., O3
  populate_gases_dry_aero(dry_aero_, buffer_);
  // cloudborne aerosol, e.g., soa_c_1
  populate_cloudborne_dry_aero(dry_aero_, buffer_);

  // FIXME: In other MAM4xx processes,
  // we are using pseudo_density instead of pseudo_density_dry to set
  // dry_atm_.p_del.
  dry_atm_.phis  = get_field_in("phis").get_view<const Real *>();
  dry_atm_.p_del = get_field_in("pseudo_density_dry").get_view<const Real **>();

  // prescribed volcanic aerosols.
  ssa_cmip6_sw_ =
      mam_coupling::view_3d("ssa_cmip6_sw", ncol_, nlev_, nswbands_);
  af_cmip6_sw_ = mam_coupling::view_3d("af_cmip6_sw", ncol_, nlev_, nswbands_);
  ext_cmip6_sw_ =
      mam_coupling::view_3d("ext_cmip6_sw", ncol_, nswbands_, nlev_);
  ext_cmip6_lw_ =
      mam_coupling::view_3d("ext_cmip6_lw_", ncol_, nlev_, nlwbands_);
  // FIXME: We need to get ssa_cmip6_sw_, af_cmip6_sw_, ext_cmip6_sw_,
  // ext_cmip6_lw_ from a netcdf file.
  // We will rely on eamxx to interpolate/map data from netcdf files.
  // The io interface in eamxx is being upgraded.
  // Thus, I will wait until changes in the eamxx's side are completed.
  Kokkos::deep_copy(ssa_cmip6_sw_, 0.0);
  Kokkos::deep_copy(af_cmip6_sw_, 0.0);
  Kokkos::deep_copy(ext_cmip6_sw_, 0.0);
  Kokkos::deep_copy(ext_cmip6_lw_, 0.0);

  const int work_len = mam4::modal_aer_opt::get_work_len_aerosol_optics();
  work_              = mam_coupling::view_2d("work", ncol_, work_len);

  // shortwave aerosol scattering asymmetry parameter [unitless]
  tau_ssa_g_sw_ =
      mam_coupling::view_3d("tau_ssa_g_sw_", ncol_, nswbands_, nlev_ + 1);
  // shortwave aerosol single-scattering albedo [unitless]
  tau_ssa_sw_ =
      mam_coupling::view_3d("tau_ssa_sw_", ncol_, nswbands_, nlev_ + 1);
  // shortwave aerosol extinction optical depth [unitless]
  tau_sw_ = mam_coupling::view_3d("tau_sw_", ncol_, nswbands_, nlev_ + 1);
  // aerosol forward scattered fraction * tau * w
  tau_f_sw_ = mam_coupling::view_3d("tau_f_sw_", ncol_, nswbands_, nlev_ + 1);

  // read table info
  {
    using namespace ShortFieldTagsNames;

    using view_1d_host = typename KT::view_1d<Real>::HostMirror;

    // Note: these functions do not set values for aerosol_optics_device_data_.
    mam4::modal_aer_opt::set_complex_views_modal_aero(
        aerosol_optics_device_data_);
    mam4::modal_aer_opt::set_aerosol_optics_data_for_modal_aero_sw_views(
        aerosol_optics_device_data_);
    mam4::modal_aer_opt::set_aerosol_optics_data_for_modal_aero_lw_views(
        aerosol_optics_device_data_);

    mam_coupling::AerosolOpticsHostData aerosol_optics_host_data;

    std::map<std::string, FieldLayout> layouts;
    std::map<std::string, view_1d_host> host_views;
    ekat::ParameterList rrtmg_params;

    mam_coupling::set_parameters_table(aerosol_optics_host_data, rrtmg_params,
                                       layouts, host_views);

    for(int imode = 0; imode < ntot_amode; imode++) {
      const auto key =
          "mam4_mode" + std::to_string(imode + 1) + "_physical_properties_file";
      const auto &fname = m_params.get<std::string>(key);
      mam_coupling::read_rrtmg_table(fname,
                                     imode,  // mode No
                                     rrtmg_params, grid_, host_views, layouts,
                                     aerosol_optics_host_data,
                                     aerosol_optics_device_data_);
    }

    std::string table_name_water =
        m_params.get<std::string>("mam4_water_refindex_file");

    // it will syn data to device.
    mam_coupling::read_water_refindex(table_name_water, grid_,
                                      aerosol_optics_device_data_.crefwlw,
                                      aerosol_optics_device_data_.crefwsw);
    //
    {
      // make a list of host views
      std::map<std::string, view_1d_host> host_views_aero;
      // defines layouts
      std::map<std::string, FieldLayout> layouts_aero;
      ekat::ParameterList params_aero;
      std::string surname_aero = "aer";
      mam_coupling::set_refindex_names(surname_aero, params_aero,
                                       host_views_aero, layouts_aero);

      constexpr int maxd_aspectype = mam4::ndrop::maxd_aspectype;
      auto specrefndxsw_host       = mam_coupling::complex_view_2d::HostMirror(
                "specrefndxsw_host", nswbands_, maxd_aspectype);

      auto specrefndxlw_host = mam_coupling::complex_view_2d::HostMirror(
          "specrefndxlw_host", nlwbands_, maxd_aspectype);

      // read physical properties data for aerosol species
      std::map<std::string, int> map_table_name_species_id;
      map_table_name_species_id["soa"]  = 4;  // soa:s-organic
      map_table_name_species_id["dust"] = 7;  // dst:dust:
      map_table_name_species_id["nacl"] = 6;  // ncl:seasalt
      map_table_name_species_id["so4"]  = 0;  // so4:sulfate
      map_table_name_species_id["pom"]  = 3;  // pom:p-organic
      map_table_name_species_id["bc"]   = 5;  // bc :black-c
      map_table_name_species_id["mom"]  = 8;  // mom:m-organic

      for(const auto &item : map_table_name_species_id) {
        const auto spec_name = item.first;
        const int species_id = item.second;
        const auto table_name =
            "mam4_" + spec_name + "_physical_properties_file";
        const auto &fname = m_params.get<std::string>(table_name);
        // read data
        // need to update table name
        params_aero.set("Filename", fname);
        AtmosphereInput refindex_aerosol(params_aero, grid_, host_views_aero,
                                         layouts_aero);
        refindex_aerosol.read_variables();
        refindex_aerosol.finalize();
        // copy data to device
        mam_coupling::set_refindex_aerosol(
            species_id, host_views_aero,
            specrefndxsw_host,  // complex refractive index for water visible
            specrefndxlw_host);
      }  // done ispec
      // reshape specrefndxsw_host and copy it to device
      mam4::modal_aer_opt::set_device_specrefindex(
          aerosol_optics_device_data_.specrefindex_sw, "short_wave",
          specrefndxsw_host);
      mam4::modal_aer_opt::set_device_specrefindex(
          aerosol_optics_device_data_.specrefindex_lw, "long_wave",
          specrefndxlw_host);
    }
  }
  // FIXME: We are hard-coding the band ordering in RRTMGP.
  // TODO: We can update optics file using the ordering below
  // (rrtmg_to_rrtmgp_swbands_).
  //  Mapping from old RRTMG sw bands to new band ordering in RRTMGP
  //  rrtmg_swband (old) = 2925, 3625, 4325, 4900, 5650, 6925, 7875, 10450,
  //  14425, 19325, 25825, 33500, 44000, 1710 ;
  // rrtmgp_swband (new) = 1710, 2925, 3625, 4325, 4900, 5650, 6925, 7875,
  // 10450, 14425, 19325, 25825, 33500, 44000 ; given the rrtmg index return the
  // rrtmgp index
  std::vector<int> temporal = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0};
  auto get_idx_rrtmgp_from_rrtmg_swbands_host =
      mam_coupling::view_int_1d::HostMirror(temporal.data(), nswbands_);
  get_idx_rrtmgp_from_rrtmg_swbands_ =
      mam_coupling::view_int_1d("rrtmg_to_rrtmgp_swbands", nswbands_);
  Kokkos::deep_copy(get_idx_rrtmgp_from_rrtmg_swbands_,
                    get_idx_rrtmgp_from_rrtmg_swbands_host);
}
void MAMOptics::run_impl(const double dt) {
  constexpr Real zero = 0.0;
  constexpr Real one  = 1.0;

  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  pre_process(wet_aero_, dry_aero_, wet_atm_, dry_atm_);
  Kokkos::fence();

  // tau_w_g : aerosol asymmetry parameter * tau * w
  const auto tau_ssa_g_sw = tau_ssa_g_sw_;
  // tau_w : aerosol single scattering albedo * tau
  const auto tau_ssa_sw = tau_ssa_sw_;
  // tau : aerosol extinction optical depth
  const auto tau_sw = tau_sw_;
  // get_field_out("aero_tau_sw_mam4").get_view<Real ***>();
  // aero_tau_lw ( or odap_aer) : absorption optical depth, per layer
  const auto aero_tau_lw = get_field_out("aero_tau_lw").get_view<Real ***>();

  const auto aero_g_sw_eamxx = get_field_out("aero_g_sw").get_view<Real ***>();

  const auto aero_ssa_sw_eamxx =
      get_field_out("aero_ssa_sw").get_view<Real ***>();
  const auto aero_tau_sw_eamxx =
      get_field_out("aero_tau_sw").get_view<Real ***>();
  // tau_w_f : aerosol forward scattered fraction * tau * w
  const auto tau_f_sw = tau_f_sw_;
  const auto aodvis   = get_field_out("aodvis").get_view<Real *>();

  // NOTE! we need a const mam_coupling::DryAtmosphere dry_atm for gpu access.
  // We cannot use member of this class inside of the parallel_for
  const mam_coupling::DryAtmosphere &dry_atm = dry_atm_;
  const auto &p_del                          = p_del_;
  const auto &ssa_cmip6_sw                   = ssa_cmip6_sw_;
  const auto &af_cmip6_sw                    = af_cmip6_sw_;
  const auto &ext_cmip6_sw                   = ext_cmip6_sw_;
  const auto &ext_cmip6_lw                   = ext_cmip6_lw_;
  const auto &work                           = work_;
  const auto &dry_aero                       = dry_aero_;
  const auto &aerosol_optics_device_data     = aerosol_optics_device_data_;
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const Int icol = team.league_rank();  // column index
        // absorption optical depth, per layer [unitless]
        const auto atm     = mam_coupling::atmosphere_for_column(dry_atm, icol);

        // FIXME: dry mass pressure interval [Pa]
        auto zi   = ekat::subview(dry_atm.z_iface, icol);
        auto pdel = ekat::subview(p_del, icol);

        auto ssa_cmip6_sw_icol = ekat::subview(ssa_cmip6_sw, icol);
        auto af_cmip6_sw_icol  = ekat::subview(af_cmip6_sw, icol);
        auto ext_cmip6_sw_icol = ekat::subview(ext_cmip6_sw, icol);

        // tau_w: aerosol single scattering albedo * tau
        auto tau_w_icol = ekat::subview(tau_ssa_sw, icol);
        // tau_w_g: aerosol assymetry
        // parameter * tau * w
        auto tau_w_g_icol = ekat::subview(tau_ssa_g_sw, icol);
        // tau_w_f: aero_tau_forward
        // forward scattered fraction * tau * w
        auto tau_w_f_icol = ekat::subview(tau_f_sw, icol);
        // tau: aerosol
        // aerosol extinction optical depth
        auto tau_icol = ekat::subview(tau_sw, icol);

        auto work_icol = ekat::subview(work, icol);

        // fetch column-specific subviews into aerosol prognostics
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);

        mam4::aer_rad_props::aer_rad_props_sw(
            team, dt, progs, atm, zi, pdel, ssa_cmip6_sw_icol, af_cmip6_sw_icol,
            ext_cmip6_sw_icol, tau_icol, tau_w_icol, tau_w_g_icol, tau_w_f_icol,
            aerosol_optics_device_data, aodvis(icol), work_icol);

      });
  Kokkos::fence();
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const Int icol = team.league_rank();  // column index
        // absorption optical depth, per layer [unitless]
        auto odap_aer_icol = ekat::subview(aero_tau_lw, icol);
        const auto atm     = mam_coupling::atmosphere_for_column(dry_atm, icol);

        // FIXME: dry mass pressure interval [Pa]
        auto zi   = ekat::subview(dry_atm.z_iface, icol);
        auto pdel = ekat::subview(p_del, icol);
        auto ext_cmip6_lw_icol = ekat::subview(ext_cmip6_lw, icol);

        // fetch column-specific subviews into aerosol prognostics
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);

        mam4::aer_rad_props::aer_rad_props_lw(
            team, dt, progs, atm, zi, pdel, ext_cmip6_lw_icol,
            aerosol_optics_device_data, odap_aer_icol);
      });
  Kokkos::fence();
  // TODO: We will need to generate optical inputs files with  band ordering
  // that is consistent with RRTMGP Optical files depend on the band ordering in
  // RRTMG. As a temporary fix, we are correcting the band ordering of mam4xx's
  // ouputs, so that they are consistent with inputs in RRTMGP. Mapping from old
  // RRTMG sw bands to new band ordering in RRTMGP than rrtmgp mam4 layout:
  // (ncols, nswlands, nlevs +1  ) rrtmgp in emaxx: (ncols, nswlands, nlevs)
  // Here, we copy data from kk=1 in mam4xx Here, we are following:
  // E3SM/components/eam/src/physics/rrtmgp
  /// cam_optics.F90
  const auto &get_idx_rrtmgp_from_rrtmg_swbands =
      get_idx_rrtmgp_from_rrtmg_swbands_;
  // postprocess output
  post_process(wet_aero_, dry_aero_, dry_atm_);
  Kokkos::fence();

  // nswbands loop is using rrtmg indexing.
  Kokkos::parallel_for(
      "copying data from mam4xx to eamxx",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0},
                                             {ncol_, nswbands_, nlev_}),
      KOKKOS_LAMBDA(const int icol, const int iswband, const int kk) {
        // Extract single scattering albedo from the product-defined fields
        if(tau_sw(icol, iswband, kk + 1) > zero) {
          aero_ssa_sw_eamxx(icol, get_idx_rrtmgp_from_rrtmg_swbands(iswband),
                            kk) =
              tau_ssa_sw(icol, iswband, kk + 1) / tau_sw(icol, iswband, kk + 1);
        } else {
          aero_ssa_sw_eamxx(icol, get_idx_rrtmgp_from_rrtmg_swbands(iswband),
                            kk) = one;
        }
        // Extract assymmetry parameter from the product-defined fields
        if(tau_ssa_sw(icol, iswband, kk + 1) > zero) {
          aero_g_sw_eamxx(icol, get_idx_rrtmgp_from_rrtmg_swbands(iswband),
                          kk) = tau_ssa_g_sw(icol, iswband, kk + 1) /
                                tau_ssa_sw(icol, iswband, kk + 1);
        } else {
          aero_g_sw_eamxx(icol, get_idx_rrtmgp_from_rrtmg_swbands(iswband),
                          kk) = zero;
        }
        // Copy cloud optical depth over directly
        aero_tau_sw_eamxx(icol, get_idx_rrtmgp_from_rrtmg_swbands(iswband),
                          kk) = tau_sw(icol, iswband, kk + 1);
      });
  Kokkos::fence();
}

void MAMOptics::finalize_impl() {}

}  // namespace scream
