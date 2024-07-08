#include <ekat/ekat_assert.hpp>
#include <physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp>

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

// for SCREAM_CIME_BUILD
#include "scream_config.h"

/*
Future work:
Write comments
write in/outs for all variables clearly
*/

namespace scream {

// =========================================================================================
MAMSrfOnlineEmiss::MAMSrfOnlineEmiss(const ekat::Comm& comm,
                                     const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// =========================================================================================
void MAMSrfOnlineEmiss::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  Units q_unit(kg / kg, "kg/kg");

  Units n_unit(1 / kg, "#/kg"); // units of number mixing ratios of tracers

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();      // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_layout_mid{{COL, LEV}, {ncol_, nlev_}};
  const FieldLayout scalar3d_layout_int{{COL, ILEV}, {ncol_, nlev_ + 1}};

  // Layout for 2D (2d horiz) variable
  const FieldLayout scalar2d_layout{{COL}, {ncol_}};

  // -------------------------------------------------------------------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // -------------------------------------------------------------------------------------------------------------------------
  add_field<Required>("T_mid", scalar3d_layout_mid, K,
                      grid_name); // temperature [K]
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa,
                      grid_name); // pressure at mid points in [Pa]
  add_field<Required>("p_int", scalar3d_layout_int, Pa,
                      grid_name); // total pressure
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,
                      grid_name); // pseudo density in [Pa]
  add_field<Required>("qv", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers"); // specific humidity
  add_field<Required>("qc", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers"); // liquid cloud water [kg/kg] wet
  add_field<Required>("qi", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers"); // ice cloud water [kg/kg] wet
  add_field<Updated>("nc", scalar3d_layout_mid, n_unit, grid_name,
                     "tracers"); // cloud liquid wet number mixing ratio
  add_field<Required>("ni", scalar3d_layout_mid, n_unit, grid_name,
                      "tracers"); // ice number mixing ratio
  add_field<Required>(
      "omega", scalar3d_layout_mid, Pa / s,
      grid_name); // Vertical pressure velocity [Pa/s] at midpoints

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);

    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name, "tracers");
    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char* int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if (strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, q_unit,
                           grid_name, "tracers");
      }
    }
  }
  // (cloud) aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char* cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);
    // printf("%s \n", int_nmr_field_name);

    add_field<Updated>(cld_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name);
    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char* cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);

      if (strlen(cld_mmr_field_name) > 0) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_layout_mid, q_unit,
                           grid_name);
      }
    }
  }

  // aerosol-related gases: mass mixing ratios
  for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char* gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_layout_mid, q_unit,
                       grid_name, "tracers");
  }
}

// =========================================================================================
// ON HOST, returns the number of bytes of device memory needed by the above
// Buffer type given the number of columns and vertical levels
size_t MAMSrfOnlineEmiss::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

// =========================================================================================
// ON HOST, initializes the Buffer type with sufficient memory to store
// intermediate (dry) quantities on the given number of columns with the given
// number of vertical levels. Returns the number of bytes allocated.
void MAMSrfOnlineEmiss::init_buffers(const ATMBufferManager& buffer_manager) {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >=
                       requested_buffer_size_in_bytes(),
                   "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(
      used_mem == requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for MAMSrfOnlineEmiss.");
}
// =========================================================================================

// inline void set_emissions_layouts(
//     const std::map<std::string, int> map_spec_id,
//     const std::string emis_type,
//     std::map<std::string, view_1d_host>& host_views,
//     mam_coupling::complex_view_2d::HostMirror&
//         specrefndxsw_host, // complex refractive index for water visible
//     mam_coupling::complex_view_2d::HostMirror& specrefndxlw_host) {

//   // names take the form "online_emis_specifier_for_SO2"
//   for (const auto& item : map_spec_id) {
//     const auto spec_name = item.first;
//     const int species_id = item.second;
//     const auto file_name = emis_type + "_emis_specifier_" + spec_name;
//     // const auto& fname = m_params.get<std::string>(file_name);
//     // update file name


//     // read data
//     AtmosphereInput srf_emissions_reader(
//         params_srf_emissions, grid_, host_views_emissions, layouts_emissions);
//     srf_emissions_reader.read_variables();
//     srf_emissions_reader.finalize();
//   } // end ispec



  // for (int i = 0; i < nswbands; i++) {
  //   specrefndxsw_host(i, species_id).real() = host_views[sw_real_name](i);
  //   specrefndxsw_host(i, species_id).imag() =
  //       haero::abs(host_views[sw_im_name](i));
  // }
  // for (int i = 0; i < nlwbands; i++) {
  //   specrefndxlw_host(i, species_id).real() = host_views[lw_real_name](i);
  //   specrefndxlw_host(i, species_id).imag() =
  //       haero::abs(host_views[lw_im_name](i));
  // }

// } // end
// =========================================================================================
void MAMSrfOnlineEmiss::initialize_impl(const RunType run_type) {
  // Gather runtime options
  //(e.g.) runtime_options.lambda_low    = m_params.get<double>("lambda_low");

  wet_atm_.qv = get_field_in("qv").get_view<const Real**>();
  wet_atm_.qc = get_field_in("qc").get_view<const Real**>();
  wet_atm_.nc = get_field_in("nc").get_view<const Real**>();
  wet_atm_.qi = get_field_in("qi").get_view<const Real**>();
  wet_atm_.ni = get_field_in("ni").get_view<const Real**>();

  dry_atm_.T_mid = get_field_in("T_mid").get_view<const Real**>();
  dry_atm_.p_mid = get_field_in("p_mid").get_view<const Real**>();
  dry_atm_.p_del = get_field_in("pseudo_density").get_view<const Real**>();
  dry_atm_.omega = get_field_in("omega").get_view<const Real**>();
  dry_atm_.qv = buffer_.qv_dry;
  dry_atm_.qc = buffer_.qc_dry;
  dry_atm_.nc = buffer_.nc_dry;
  dry_atm_.qi = buffer_.qi_dry;
  dry_atm_.ni = buffer_.ni_dry;

  // NOTE: these are taken as arguments to srf_emissions_inti()
  // and then passed to trcdata_init()
  // rmv_file = false;
  // emis_cycle_yr
  // emis_fixed_ymd
  // emis_fixed_tod
  // emis_type

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    wet_aero_.int_aero_nmr[m] =
        get_field_out(int_nmr_field_name).get_view<Real**>();
    dry_aero_.int_aero_nmr[m] = buffer_.dry_int_aero_nmr[m];

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    const char* cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);
    wet_aero_.cld_aero_nmr[m] =
        get_field_out(cld_nmr_field_name).get_view<Real**>();
    dry_aero_.cld_aero_nmr[m] = buffer_.dry_cld_aero_nmr[m];

    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char* int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);
      if (strlen(int_mmr_field_name) > 0) {
        wet_aero_.int_aero_mmr[m][a] =
            get_field_out(int_mmr_field_name).get_view<Real**>();
        dry_aero_.int_aero_mmr[m][a] = buffer_.dry_int_aero_mmr[m][a];
      }

      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const char* cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if (strlen(cld_mmr_field_name) > 0) {
        wet_aero_.cld_aero_mmr[m][a] =
            get_field_out(cld_mmr_field_name).get_view<Real**>();
        dry_aero_.cld_aero_mmr[m][a] = buffer_.dry_cld_aero_mmr[m][a];
      }
    }
  }
  for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char* gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] = get_field_out(gas_mmr_field_name).get_view<Real**>();
    dry_aero_.gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

  // set up our preprocess functor
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                         dry_aero_);

  // // read data from file
  {
  //   // using namespace ShortFieldTagsNames;

    using view_1d_host = typename KT::view_1d<Real>::HostMirror;

    mam_coupling::AerosolSurfaceEmissionsHostData srf_emissions_host_data;

  //   {
      // make a list of host views
      std::map<std::string, view_1d_host> host_views_srf_emissions;
      // defines layouts
      std::map<std::string, FieldLayout> layouts_srf_emissions;
      ekat::ParameterList params_srf_emissions;
      std::string prefix_srf_emissions = "srf_emis_specifier_for_";

  //     // constexpr int maxd_aspectype = mam4::ndrop::maxd_aspectype;
  //     // auto specrefndxsw_host = mam_coupling::complex_view_2d::HostMirror(
  //     //     "specrefndxsw_host", nswbands_, maxd_aspectype);

  //     // auto specrefndxlw_host = mam_coupling::complex_view_2d::HostMirror(
  //     //     "specrefndxlw_host", nlwbands_, maxd_aspectype);

      std::map<std::string, int> map_srf_emiss_name_species_id;
      map_srf_emiss_name_species_id["DMS"] = 0;
      map_srf_emiss_name_species_id["SO2"] = 1;
      map_srf_emiss_name_species_id["bc_a4"] = 2;
      map_srf_emiss_name_species_id["num_a1"] = 3;
      map_srf_emiss_name_species_id["num_a2"] = 4;
      map_srf_emiss_name_species_id["num_a4"] = 5;
      map_srf_emiss_name_species_id["pom_a4"] = 6;
      map_srf_emiss_name_species_id["so4_a1"] = 7;
      map_srf_emiss_name_species_id["so4_a2"] = 8;

      std::map<std::string, int> map_online_emiss_name_species_id;
      map_online_emiss_name_species_id["SO2"] = 0;
      map_online_emiss_name_species_id["SOAG"] = 1;
      map_online_emiss_name_species_id["bc_a4"] = 2;
      map_online_emiss_name_species_id["num_a1"] = 3;
      map_online_emiss_name_species_id["num_a2"] = 4;
      map_online_emiss_name_species_id["num_a4"] = 5;
      map_online_emiss_name_species_id["pom_a4"] = 6;
      map_online_emiss_name_species_id["so4_a1"] = 7;
      map_online_emiss_name_species_id["so4_a2"] = 8;

      // To create the input object, we need to set:
      // 1) names (do during read loop)
      // 2) params
      // 3) host views
      // 4) layouts
      // set_emissions_names(surname_emissions, params_srf_emissions,
      // host_views_emissions, layouts_emissions);

      // inline void set_emissions_names(
      //     const std::map<std::string, int> map_spec_id,
      //     const std::string emis_type,
      //     std::map<std::string, view_1d_host>& host_views,
      //     mam_coupling::complex_view_2d::HostMirror&
      //         specrefndxsw_host, // complex refractive index for water visible
      //     mam_coupling::complex_view_2d::HostMirror& specrefndxlw_host) {

      //   // names take the form "online_emis_specifier_for_SO2"
      //   for (const auto& item : map_spec_id) {
      //     const auto spec_name = item.first;
      //     const int species_id = item.second;
      //     const auto file_name = emis_type + "_emis_specifier_" + spec_name;
          // const auto& fname = m_params.get<std::string>(file_name);
          // update file name


          // // read data
          // AtmosphereInput srf_emissions_reader(
          //     params_srf_emissions, grid_, host_views_emissions, layouts_emissions);
          // srf_emissions_reader.read_variables();
          // srf_emissions_reader.finalize();
        // } // end ispec



        // for (int i = 0; i < nswbands; i++) {
        //   specrefndxsw_host(i, species_id).real() = host_views[sw_real_name](i);
        //   specrefndxsw_host(i, species_id).imag() =
        //       haero::abs(host_views[sw_im_name](i));
        // }
        // for (int i = 0; i < nlwbands; i++) {
        //   specrefndxlw_host(i, species_id).real() = host_views[lw_real_name](i);
        //   specrefndxlw_host(i, species_id).imag() =
        //       haero::abs(host_views[lw_im_name](i));
        // }

      // } // end

      // ===============
      // Names
      // ===============

      // names take the form <type>_emis_specifier_for_<spec>
      // set_emissions_names("srf", params_srf_emissions, host_views_srf_emissions, layouts_srf_emissions);
      // set_emissions_names("online", params_srf_emissions, host_views_online_emissions, layouts_online_emissions);

      // ===============
      // Params
      // ===============
      using strvec_t = std::vector<std::string>;
      params_srf_emissions.set("Skip_Grid_Checks", true);
      // params_srf_emissions.set<strvec_t>("Field Names", {refindex_real_sw, refindex_im_sw,
      //                                      refindex_real_lw, refindex_im_lw});
      // params_srf_emissions.set("Filename", fname);

      // ===============
      // Host Views
      // ===============
      using view_1d_host = typename KT::view_1d<Real>::HostMirror;
      std::map<std::string, view_1d_host> host_views;

      // host_views[refindex_real_sw] = view_1d_host(refindex_real_sw,
      // nswbands); host_views[refindex_im_sw]   = view_1d_host(refindex_im_sw,
      // nswbands); host_views[refindex_real_lw] =
      // view_1d_host(refindex_real_lw, nlwbands); host_views[refindex_im_lw] =
      // view_1d_host(refindex_im_lw, nlwbands);

      // ===============
      // Layouts
      // ===============
      std::map<std::string, FieldLayout> layouts;
      // FieldLayout scalar_refindex_sw_layout{{SWBND}, {nswbands}};
      // FieldLayout scalar_refindex_lw_layout{{LWBND}, {nlwbands}};

      // layouts.emplace(refindex_real_sw, scalar_refindex_sw_layout);
      // layouts.emplace(refindex_im_sw, scalar_refindex_sw_layout);
      // layouts.emplace(refindex_real_lw, scalar_refindex_lw_layout);
      // layouts.emplace(refindex_im_lw, scalar_refindex_lw_layout);

      // reshape specrefndxsw_host and copy it to device
      // mam4::modal_aer_opt::set_device_specrefindex(
      //     aerosol_optics_device_data_.specrefindex_sw, "short_wave",
      //     specrefndxsw_host);
      // mam4::modal_aer_opt::set_device_specrefindex(
      //     aerosol_optics_device_data_.specrefindex_lw, "long_wave",
      //     specrefndxlw_host);
  }
}

// =============================================================================
void MAMSrfOnlineEmiss::run_impl(const double dt) {

  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  /* Rough notes:

  Here we should implement or port the chem_emissions subroutine in
  chemistry.F90. Basically call two subroutines, aero_model_emissions and
  set_srf_emissions.

  Here is the code:

    ! initialize chemistry constituent surface fluxes to zero
    do m = 2,pcnst
       n = map2chm(m)
       if (n>0) cam_in%cflx(:,m) = 0._r8
    enddo

    ! aerosol emissions ...
    call aero_model_emissions( state, & ! in
                               cam_in ) ! out

    ! prescribed emissions from file ...

    !-----------------------------------------------------------------------
    !        ... Set surface emissions
    !-----------------------------------------------------------------------
    call set_srf_emissions( lchnk, ncol, sflx(:,:) )

    do m = 1,pcnst
       n = map2chm(m)
       if ( n /= h2o_ndx .and. n > 0 ) then
          cam_in%cflx(:ncol,m) = cam_in%cflx(:ncol,m) + sflx(:ncol,n)
          call outfld( sflxnam(m), cam_in%cflx(:ncol,m), ncol,lchnk )
       endif
    enddo


  */

  std::cout << "End of surface emissions run" << std::endl;
}

// =============================================================================
} // namespace scream
