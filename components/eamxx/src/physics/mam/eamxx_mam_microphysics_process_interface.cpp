#include <mam4xx/mam4.hpp>
#include <physics/mam/eamxx_mam_microphysics_process_interface.hpp>
// impl namespace for some driver level functions for microphysics

#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"
#include "readfiles/find_season_index_utils.hpp"
#include "readfiles/photo_table_utils.cpp"
#include "physics/mam/readfiles/vertical_remapper_mam4.hpp"
#include "share/util/eamxx_data_interpolation.hpp"

namespace scream {

MAMMicrophysics::MAMMicrophysics(const ekat::Comm &comm,
                                 const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params), aero_config_() {
  config_.amicphys.do_cond   = m_params.get<bool>("mam4_do_cond");
  config_.amicphys.do_rename = m_params.get<bool>("mam4_do_rename");
  config_.amicphys.do_newnuc = m_params.get<bool>("mam4_do_newnuc");
  config_.amicphys.do_coag   = m_params.get<bool>("mam4_do_coag");
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);

  // these parameters guide the coupling between parameterizations
  // NOTE: mam4xx was ported with these parameters fixed, so it's probably not
  // NOTE: safe to change these without code modifications.

  // controls treatment of h2so4 condensation in mam_gasaerexch_1subarea
  //    1 = sequential   calc. of gas-chem prod then condensation loss
  //    2 = simultaneous calc. of gas-chem prod and  condensation loss
  config_.amicphys.gaexch_h2so4_uptake_optaa = 2;

  // controls treatment of h2so4 concentrationin mam_newnuc_1subarea
  //    1 = use average value calculated in standard cam5.2.10 and earlier
  //    2 = use average value calculated in mam_gasaerexch_1subarea
  //   11 = use average of initial and final values from mam_gasaerexch_1subarea
  //   12 = use final value from mam_gasaerexch_1subarea
  config_.amicphys.newnuc_h2so4_conc_optaa = 2;

  // LINOZ namelist parameters
  // Compute LINOZ only for prognostic Ozone (i.e. !use_prescribed_ozone_)
  config_.linoz.compute = !use_prescribed_ozone_;

  if (config_.linoz.compute) {
    config_.linoz.o3_lbl = m_params.get<int>("mam4_o3_lbl");
    config_.linoz.o3_tau = m_params.get<double>("mam4_o3_tau");
    config_.linoz.o3_sfc = m_params.get<double>("mam4_o3_sfc");
    config_.linoz.psc_T  = m_params.get<double>("mam4_psc_T");
  }

}
// ================================================================
//  SET_GRIDS
// ================================================================
void MAMMicrophysics::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  // set grid for all the inputs and outputs
  // use physics grid
  grid_                 = grids_manager->get_grid("physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // number of levels per column

  // create our photolysis rate calculation table
  const std::string rsf_file = m_params.get<std::string>("mam4_rsf_file");
  const std::string xs_long_file =
      m_params.get<std::string>("mam4_xs_long_file");

  photo_table_ = impl::read_photo_table(rsf_file, xs_long_file);
  // NOTE: we need photo_table_ before getting len_temporary_views_.
  len_temporary_views_ = get_len_temporary_views();
  buffer_.set_len_temporary_views(len_temporary_views_);
  buffer_.set_num_scratch(num_2d_scratch_);

  // get column geometry and locations
  col_latitudes_ = grid_->get_geometry_data("lat").get_view<const Real *>();

  // define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (2d horiz) variable
  const FieldLayout scalar2d = grid_->get_2d_scalar_layout();

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);
  const FieldLayout scalar3d_int = grid_->get_3d_scalar_layout(false);

  // For U and V components of wind
  const FieldLayout vector3d = grid_->get_3d_vector_layout(true, 2);

  using namespace ekat::units;
  // --------------------------------------------------------------------------
  // These variables are "Required" or pure inputs for the process
  // --------------------------------------------------------------------------

  add_tracers_wet_atm();
  add_fields_dry_atm();

  // cloud liquid number mixing ratio [1/kg]
  auto n_unit           = 1 / kg;   // units of number mixing ratios of tracers
  add_tracer<Required>("nc", grid_, n_unit);

  constexpr auto m2 = pow(m, 2);
  constexpr auto s2 = pow(s, 2);

  // Surface geopotential [m2/s2]
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  // Surface pressure [Pa]
  add_field<Required>("ps", scalar2d, Pa, grid_name);

  // U and V components of the wind[m/s]
  add_field<Required>("horiz_winds", vector3d, m / s, grid_name);

  //----------- Variables from microphysics scheme -------------
  constexpr auto nondim = ekat::units::Units::nondimensional();
  // Total cloud fraction [fraction]
  add_field<Required>("cldfrac_liq", scalar3d_mid, nondim, grid_name);

  // Evaporation from stratiform rain [kg/kg/s]
  add_field<Required>("nevapr", scalar3d_mid, kg / kg / s, grid_name);

  // Stratiform rain production rate [kg/kg/s]
  add_field<Required>("precip_total_tend", scalar3d_mid, kg / kg / s,
                      grid_name);

  // precipitation liquid mass [kg/m2]
  add_field<Required>("precip_liq_surf_mass", scalar2d, kg / m2, grid_name);

  // precipitation ice mass [kg/m2]
  add_field<Required>("precip_ice_surf_mass", scalar2d, kg / m2, grid_name);

  //----------- Variables from other mam4xx processes ------------
  // Number of modes
  constexpr int nmodes = mam4::AeroConfig::num_modes();

  // layout for 3D (ncol, nmodes, nlevs)
  FieldLayout vector3d_mid_nmodes = grid_->get_3d_vector_layout(
      true, nmodes, mam_coupling::num_modes_tag_name());

  // Geometric mean dry diameter for number distribution [m]
  add_field<Required>("dgnum", vector3d_mid_nmodes, m, grid_name);
  // Geometric mean wet diameter for number distribution [m]
  add_field<Required>("dgnumwet", vector3d_mid_nmodes, m, grid_name);

  constexpr auto m3 = pow(m, 3);
  // Wet density of interstitial aerosol [kg/m3]
  add_field<Required>("wetdens", vector3d_mid_nmodes, kg / m3, grid_name);

  // For fractional land use
  const FieldLayout vector2d_class =
      grid_->get_2d_vector_layout(mam4::mo_drydep::n_land_type, "class");

  // Fractional land use [fraction]
  add_field<Required>("fraction_landuse", vector2d_class, nondim, grid_name);

  //----------- Variables from the coupler ---------
  // surface albedo shortwave, direct
  add_field<Required>("sfc_alb_dir_vis", scalar2d, nondim, grid_name);

  // Surface temperature[K]
  add_field<Required>("surf_radiative_T", scalar2d, K, grid_name);

  // snow depth land [m]
  add_field<Required>("snow_depth_land", scalar2d, m, grid_name);

  //----------- Variables from the RRTMGP radiation ---------
  // Downwelling solar flux at the surface [w/m2]
  add_field<Required>("SW_flux_dn", scalar3d_int, W / m2, grid_name);

  // ---------------------------------------------------------------------
  // These variables are "updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  // add tracers, e.g., num_a1, soa_a1
  add_tracers_interstitial_aerosol();
  // add tracer gases, e.g., O3
  add_tracers_gases();
  // add fields e.g., num_c1, soa_c1
  add_fields_cloudborne_aerosol();
  //----------- Updated variables from other mam4xx processes ------------
  // layout for Constituent fluxes
  FieldLayout vector2d_pcnst =
      grid_->get_2d_vector_layout(mam4::pcnst, "num_phys_constituents");
  // Constituent fluxes of species in [kg/m2/s]
  add_field<Updated>("constituent_fluxes", vector2d_pcnst, kg / m2 / s,
                     grid_name);

  // ---------------------------------------------------------------------
  // These variables are "computed" or outputs for the process
  // ---------------------------------------------------------------------
  // Number of externally forced chemical species
  constexpr int extcnt = mam4::gas_chemistry::extcnt;

  FieldLayout vector3d_extcnt = grid_->get_3d_vector_layout(true, extcnt, "ext_cnt");

  // Register computed fields for external forcing
  // - extfrc: 3D instantaneous forcing rate [kg/m³/s]
  add_field<Computed>("mam4_external_forcing", vector3d_extcnt, kg / m3 / s, grid_name);

  // Diagnostic fields for aerosol microphysics

  //Flag to indicate if we want to compute extra diagnostics
  extra_mam4_aero_microphys_diags_ = m_params.get<bool>("extra_mam4_aero_microphys_diags", false);
  if (extra_mam4_aero_microphys_diags_) {
    const FieldLayout vector3d_num_gas_aerosol_constituents =
        grid_->get_3d_vector_layout(true, mam_coupling::gas_pcnst(), "num_gas_aerosol_constituents");

    const FieldLayout vector2d_nmodes =
      grid_->get_2d_vector_layout(nmodes, "nmodes");

    // Diagnostics: tendencies due to gas phase chemistry [mixed units: kg/kg/s or #/kg/s]
    add_field<Computed>("mam4_microphysics_tendency_gas_phase_chemistry", vector3d_num_gas_aerosol_constituents, nondim, grid_name);

    // Diagnostics: tendencies due to aqueous chemistry [mixed units: kg/kg/s or #/kg/s]
    add_field<Computed>("mam4_microphysics_tendency_aqueous_chemistry", vector3d_num_gas_aerosol_constituents, nondim, grid_name);

    // Diagnostics: SO4 in-cloud tendencies [mixed units: kg/kg/s or #/kg/s]
    add_field<Computed>("mam4_microphysics_tendency_aqso4", vector3d_mid_nmodes, nondim, grid_name);

    // Diagnostics: H2SO4 in-cloud tendencies [mixed units: kg/kg/s or #/kg/s]
    add_field<Computed>("mam4_microphysics_tendency_aqh2so4", vector3d_mid_nmodes, nondim, grid_name);

    // Register computed diagnostic fields
    //(NOTE: dqdt_so4_aqueous_chemistry is the vertically reduced field of "mam4_microphysics_tendency_aqso4")
    add_field<Computed>("dqdt_so4_aqueous_chemistry", vector2d_nmodes, kg/m2/s,  grid_name);
    //(NOTE: dqdt_h2so4_uptake is the vertically reduced field of "mam4_microphysics_tendency_aqh2so4")
    add_field<Computed>("dqdt_h2so4_uptake", vector2d_nmodes, kg/m2/s,  grid_name);

    // Diagnostics: tendencies due to aerosol microphysics (gas aerosol exchange) [mixed units: mol/mol/s or #/mol/s]
    add_field<Computed>("mam4_microphysics_tendency_condensation", vector3d_num_gas_aerosol_constituents, nondim, grid_name);
    add_field<Computed>("mam4_microphysics_tendency_renaming", vector3d_num_gas_aerosol_constituents, nondim, grid_name);
    add_field<Computed>("mam4_microphysics_tendency_nucleation", vector3d_num_gas_aerosol_constituents, nondim, grid_name);
    add_field<Computed>("mam4_microphysics_tendency_coagulation", vector3d_num_gas_aerosol_constituents, nondim, grid_name);
    add_field<Computed>("mam4_microphysics_tendency_renaming_cloud_borne", vector3d_num_gas_aerosol_constituents, nondim, grid_name);
  }

<<<<<<< HEAD
  // Creating a Linoz reader and setting Linoz parameters involves reading data
  // from a file and configuring the necessary parameters for the Linoz model.
#ifdef USE_OLD_LINOZ_FILE_READ
  if (config_.linoz.compute) {

=======
>>>>>>> b16a6664e4 (EAMxx: Fixing variable names and add comments.)
  // NOTE: order matches mam4xx:
  // names of variables in oxid file.
  var_names_oxi_={"O3", "OH", "NO3", "HO2"};
  // names of variables for linoz
  var_names_linoz_ = {
        "o3_clim",  "o3col_clim", "t_clim",      "PmL_clim",
        "dPmL_dO3", "dPmL_dT",    "dPmL_dO3col", "cariolle_pscs"};
<<<<<<< HEAD
    std::cout << "Using former reader..." << "\n";
=======

#ifdef USE_OLD_LINOZ_FILE_READ
  {
>>>>>>> b16a6664e4 (EAMxx: Fixing variable names and add comments.)
    linoz_file_name_ = m_params.get<std::string>("mam4_linoz_file_name");
    const std::string linoz_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");
    // in format YYYYMMDD
    const int linoz_cyclical_ymd = m_params.get<int>("mam4_linoz_ymd");
    scream::mam_coupling::setup_tracer_data(linoz_data_, linoz_file_name_,
                                            linoz_cyclical_ymd);
    LinozHorizInterp_ = scream::mam_coupling::create_horiz_remapper(
        grid_, linoz_file_name_, linoz_map_file, var_names_linoz_, linoz_data_);
    LinozDataReader_ = scream::mam_coupling::create_tracer_data_reader(
        LinozHorizInterp_, linoz_file_name_);

    // linoz reader
    const auto io_grid_linoz = LinozHorizInterp_->get_tgt_grid();
    const int num_cols_io_linoz =
        io_grid_linoz->get_num_local_dofs();  // Number of columns on this rank
    const int num_levs_io_linoz =
        io_grid_linoz
            ->get_num_vertical_levels();  // Number of levels per column
    const int nvars = int(var_names_linoz_.size());
    linoz_data_.init(num_cols_io_linoz, num_levs_io_linoz, nvars);
    linoz_data_.allocate_temporary_views();
  }  // LINOZ reader

  {
    oxid_file_name_ = m_params.get<std::string>("mam4_oxid_file_name");
    const std::string oxid_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");

    // in format YYYYMMDD
    const int oxid_ymd = m_params.get<int>("mam4_oxid_ymd");
    scream::mam_coupling::setup_tracer_data(tracer_data_, oxid_file_name_,
                                            oxid_ymd);
    TracerHorizInterp_ = scream::mam_coupling::create_horiz_remapper(
        grid_, oxid_file_name_, oxid_map_file, var_names_oxi_, tracer_data_);
    TracerDataReader_ = scream::mam_coupling::create_tracer_data_reader(
        TracerHorizInterp_, oxid_file_name_);

    const int nvars    = int(var_names_oxi_.size());
    const auto io_grid = TracerHorizInterp_->get_tgt_grid();
    const int num_cols_io =
        io_grid->get_num_local_dofs();  // Number of columns on this rank
    const int num_levs_io =
        io_grid->get_num_vertical_levels();  // Number of levels per column
    tracer_data_.init(num_cols_io, num_levs_io, nvars);
    tracer_data_.allocate_temporary_views();

    for(int ivar = 0; ivar < nvars; ++ivar) {
      cnst_offline_[ivar] = view_2d("cnst_offline_", ncol_, nlev_);
    }
  }  // oxid file reader

#else
  // The DataInterpolation class uses Field. We save these fields in FM.
  for(const auto &field_name : var_names_linoz_) {
      add_field<Computed>(field_name, scalar3d_mid, nondim, grid_name);
  }
  for(const auto &field_name : var_names_oxi_) {
      add_field<Computed>(field_name, scalar3d_mid, nondim, grid_name);
  }
#endif
  // list of species for elevated emissiones.
  extfrc_lst_ = {"so2",    "so4_a1", "so4_a2", "pom_a4", "bc_a4",
                   "num_a1", "num_a2", "num_a4", "soag"};
  // sectors or variables for each species.
  elevated_emis_var_names_["so2"]    = {"BB", "ENE_ELEV", "IND_ELEV",
                                          "contvolc"};
  elevated_emis_var_names_["so4_a1"] = {"BB", "ENE_ELEV", "IND_ELEV",
                                          "contvolc"};
  elevated_emis_var_names_["so4_a2"] = {"contvolc"};
  elevated_emis_var_names_["pom_a4"] = {"BB"};
  elevated_emis_var_names_["bc_a4"]  = {"BB"};
  elevated_emis_var_names_["num_a1"] = {
        "num_a1_SO4_ELEV_BB", "num_a1_SO4_ELEV_ENE", "num_a1_SO4_ELEV_IND",
        "num_a1_SO4_ELEV_contvolc"};
  elevated_emis_var_names_["num_a2"] = {"num_a2_SO4_ELEV_contvolc"};
  // num_a4
  // FIXME: why the sectors in this files are num_a1;
  //  I guess this should be num_a4? Is this a bug in the orginal nc files?
  elevated_emis_var_names_["num_a4"] = {"num_a1_BC_ELEV_BB",
                                          "num_a1_POM_ELEV_BB"};
  elevated_emis_var_names_["soag"] = {"SOAbb_src", "SOAbg_src", "SOAff_src"};

#ifdef USE_OLD_VERTICAL_FILE_READ
  {
    const std::string extfrc_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");
    // NOTE: order of forcing species is important.
    // extfrc_lst(:  9) = {'SO2             ','so4_a1          ','so4_a2
    // ','pom_a4          ','bc_a4           ', 'num_a1          ','num_a2
    // ','num_a4          ','SOAG            ' }
    // This order corresponds to files in namelist e3smv2
    for(const auto &var_name : extfrc_lst_) {
      std::string item_name = "mam4_" + var_name + "_elevated_emiss_file_name";
      const auto file_name  = m_params.get<std::string>(item_name);
      elevated_emis_file_name_[var_name] = file_name;
    }
    int elevated_emiss_cyclical_ymd = m_params.get<int>("elevated_emiss_ymd");

    for(const auto &var_name : extfrc_lst_) {
      const auto file_name = elevated_emis_file_name_[var_name];
      const auto var_names = elevated_emis_var_names_[var_name];

      scream::mam_coupling::TracerData data_tracer;
      scream::mam_coupling::setup_tracer_data(data_tracer, file_name,
                                              elevated_emiss_cyclical_ymd);
      auto hor_rem = scream::mam_coupling::create_horiz_remapper(
          grid_, file_name, extfrc_map_file, var_names, data_tracer);

      auto file_reader = scream::mam_coupling::create_tracer_data_reader(
          hor_rem, file_name, data_tracer.file_type);
      ElevatedEmissionsHorizInterp_.push_back(hor_rem);
      ElevatedEmissionsDataReader_.push_back(file_reader);
      elevated_emis_data_.push_back(data_tracer);
    }  // var_name elevated emissions
    int i               = 0;
    for(const auto &var_name : extfrc_lst_) {
      const auto file_name = elevated_emis_file_name_[var_name];
      const auto var_names = elevated_emis_var_names_[var_name];
      const int nvars      = static_cast<int>(var_names.size());

      forcings_[i].nsectors = nvars;
      // I am assuming the order of species in extfrc_lst_.
      // Indexing in mam4xx is fortran.
      forcings_[i].frc_ndx = i + 1;
      const auto io_grid_emis =
          ElevatedEmissionsHorizInterp_[i]->get_tgt_grid();
      const int num_cols_io_emis =
          io_grid_emis->get_num_local_dofs();  // Number of columns on this rank
      const int num_levs_io_emis =
          io_grid_emis
              ->get_num_vertical_levels();  // Number of levels per column
      elevated_emis_data_[i].init(num_cols_io_emis, num_levs_io_emis, nvars);
      elevated_emis_data_[i].allocate_temporary_views();
      forcings_[i].file_alt_data = elevated_emis_data_[i].has_altitude_;
      EKAT_REQUIRE_MSG(
        nvars <= int(mam_coupling::MAX_SECTION_NUM_FORCING),
        "Error! Number of sections is bigger than "
        "MAX_SECTION_NUM_FORCING. Increase the "
        "MAX_SECTION_NUM_FORCING in tracer_reader_utils.hpp \n");
      for(int isp = 0; isp < nvars; ++isp) {
        forcings_[i].fields[isp] =
            view_2d("elevated_emis_output_", ncol_, nlev_);
      }
      ++i;
    }  // end i

  }  // Tracer external forcing data
#else
  // Fields for elevated emissions.
  for(const auto &pair : elevated_emis_var_names_) {
    const auto &var_name =pair.first;
    for(const auto &field_name : pair.second) {
      add_field<Computed>(field_name+"_"+var_name, scalar3d_mid, nondim, grid_name);
    }
  }

#endif

  {
    const std::string season_wes_file =
        m_params.get<std::string>("mam4_season_wes_file");
    const auto &clat = col_latitudes_;
    mam_coupling::find_season_index_reader(season_wes_file, clat,
                                           index_season_lai_);
  }

}  // set_grids

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by
// the above. Buffer type given the number of columns and vertical
// levels
size_t MAMMicrophysics::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_, num_2d_scratch_,
                                   len_temporary_views_);
}

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to
// store intermediate (dry) quantities on the given number of
// columns with the given number of vertical levels. Returns the
// number of bytes allocated.

void MAMMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) {
  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMMicrophysics."
                   " Used memory: "
                       << std::to_string(used_mem)
                       << "."
                          " Requested memory: "
                       << std::to_string(requested_buffer_size_in_bytes())
                       << ". \n");
}
int MAMMicrophysics::get_len_temporary_views() {
  const int photo_table_len = get_photo_table_work_len(photo_table_);
  const int sethet_work_len = mam4::mo_sethet::get_total_work_len_sethet();
  constexpr int extcnt      = mam4::gas_chemistry::extcnt;
  int work_len              = 0;
  // work_photo_table_
  work_len += ncol_ * photo_table_len;
  // work_set_het_
  work_len += ncol_ * sethet_work_len;
  // photo_rates_
  work_len += ncol_ * nlev_ * mam4::mo_photo::phtcnt;
  // invariants_
  work_len += ncol_ * nlev_ * mam4::gas_chemistry::nfs;
  // extfrc_
  work_len += ncol_ * nlev_ * extcnt;
  return work_len;
}
void MAMMicrophysics::init_temporary_views() {
  const int photo_table_len = get_photo_table_work_len(photo_table_);
  const int sethet_work_len = mam4::mo_sethet::get_total_work_len_sethet();
  constexpr int extcnt      = mam4::gas_chemistry::extcnt;
  auto work_ptr             = (Real *)buffer_.temporary_views.data();

  work_photo_table_ = view_2d(work_ptr, ncol_, photo_table_len);
  work_ptr += ncol_ * photo_table_len;
  work_set_het_ = view_2d(work_ptr, ncol_, sethet_work_len);
  work_ptr += ncol_ * sethet_work_len;
  // here's where we store per-column photolysis rates
  photo_rates_ = view_3d(work_ptr, ncol_, nlev_, mam4::mo_photo::phtcnt);
  work_ptr += ncol_ * nlev_ * mam4::mo_photo::phtcnt;
  invariants_ = view_3d(work_ptr, ncol_, nlev_, mam4::gas_chemistry::nfs);
  work_ptr += ncol_ * nlev_ * mam4::gas_chemistry::nfs;
  extfrc_ = view_3d(work_ptr, ncol_, nlev_, extcnt);
  work_ptr += ncol_ * nlev_ * extcnt;

  // Error check
  // NOTE: workspace_provided can be larger than workspace_used, but let's try
  // to use the minimum amount of memory
  const int workspace_used     = work_ptr - buffer_.temporary_views.data();
  const int workspace_provided = buffer_.temporary_views.extent(0);
  EKAT_REQUIRE_MSG(workspace_used == workspace_provided,
                   "Error: workspace_used (" + std::to_string(workspace_used) +
                       ") and workspace_provided (" +
                       std::to_string(workspace_provided) +
                       ") should be equal. \n");
}
#ifndef USE_OLD_LINOZ_FILE_READ
// set DataInterpolation object for oxid reader.
void MAMMicrophysics::set_oxid_reader()
{

  auto pmid = get_field_in("p_mid");
  // Oxid fields read initialization
  const auto oxid_file_name = m_params.get<std::string>("mam4_oxid_file_name");
  const std::string oxid_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");
  // get fields from FM.
  std::vector<Field> oxid_fields;
  for(const auto &field_name : var_names_oxi_) {
      oxid_fields.push_back(get_field_out(field_name));
  }

  const int oxid_ymd = m_params.get<int>("mam4_oxid_ymd");
  util::TimeStamp ref_ts_oxid= mam_coupling::convert_date(oxid_ymd);
  data_interp_oxid_ = std::make_shared<DataInterpolation>(grid_,oxid_fields);
  data_interp_oxid_->setup_time_database ({oxid_file_name},util::TimeLine::YearlyPeriodic, ref_ts_oxid);
  data_interp_oxid_->create_horiz_remappers (oxid_map_file=="none" ? "" : oxid_map_file);

  DataInterpolation::VertRemapData remap_data_oxid;
  remap_data_oxid.vr_type = DataInterpolation::Dynamic3DRef;
  remap_data_oxid.pname = "PS";
  remap_data_oxid.pmid = pmid;
  //FIXME: make it a namelist parameter.
  // Dynamic3DRef can also be employed instead of mam4xx routine.
  bool mam4_use_mam4xx_oxi_vert_remap=true;
  if (mam4_use_mam4xx_oxi_vert_remap){
    // We are using a custom remapper that invokes the MAM4XX routine
    // for vertical interpolation.
    // The type used is VertRemapType::MAM4_PSRef.
    auto grid_after_hremap = data_interp_oxid_->get_grid_after_hremap();
    auto vertical_remapper= std::make_shared<VerticalRemapperMAM4>(grid_after_hremap, grid_,
    VerticalRemapperMAM4::VertRemapType::MAM4_PSRef);
    remap_data_oxid.custom_remapper=vertical_remapper;
  }
  data_interp_oxid_->create_vert_remapper (remap_data);
  data_interp_oxid_->init_data_interval (start_of_step_ts());
}
// set DataInterpolation object for linoz reader.
void MAMMicrophysics::set_linoz_reader(){
  auto pmid = get_field_in("p_mid");
  const int linoz_cyclical_ymd = m_params.get<int>("mam4_linoz_ymd");
  // in format YYYYMMDD
  util::TimeStamp ref_ts_linoz = mam_coupling::convert_date(linoz_cyclical_ymd);
  const auto m_linoz_file_name = m_params.get<std::string>("mam4_linoz_file_name");
  const std::string linoz_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");
  std::vector<Field> linoz_fields;
  for(const auto &field_name : var_names_linoz_) {
      linoz_fields.push_back(get_field_out(field_name));
  }

  data_interp_linoz_ = std::make_shared<DataInterpolation>(grid_,linoz_fields);
  data_interp_linoz_->setup_time_database ({m_linoz_file_name},util::TimeLine::YearlyPeriodic, ref_ts_linoz);
  data_interp_linoz_->create_horiz_remappers (linoz_map_file=="none" ? "" : linoz_map_file);

  DataInterpolation::VertRemapData remap_data_linoz;
  remap_data_linoz.vr_type = DataInterpolation::Static1D;
  // lev is the name of variables for vertical interpolation
  remap_data_linoz.pname = "lev";
  remap_data_linoz.pmid = pmid;
  // Static1D can also be employed instead of mam4xx routine.
  bool mam4_use_mam4xx_linoz_vert_remap=true;
  if (mam4_use_mam4xx_linoz_vert_remap){
    // We are using a custom remapper that invokes the MAM4XX routine
    // for vertical interpolation.
    // The type used is VertRemapType::MAM4_ZONAL.
    auto grid_after_hremap_linoz = data_interp_linoz_->get_grid_after_hremap();
    auto vertical_remapper_linoz = std::make_shared<VerticalRemapperMAM4>(grid_after_hremap_linoz, grid_,
    VerticalRemapperMAM4::VertRemapType::MAM4_ZONAL);
    remap_data_linoz.custom_remapper=vertical_remapper_linoz;
  }
  data_interp_linoz_->create_vert_remapper (remap_data_linoz);
  data_interp_linoz_->init_data_interval (start_of_step_ts());
}
#endif

#ifndef USE_OLD_VERTICAL_FILE_READ
// set DataInterpolation object for elevated emissions reader.
void MAMMicrophysics::set_elevated_emissions_reader()
{
  const auto z_iface = get_field_out("z_iface");
  const std::string extfrc_map_file =
        m_params.get<std::string>("aero_microphys_remap_file", "");
  int elevated_emiss_cyclical_ymd = m_params.get<int>("elevated_emiss_ymd");
  for(const auto &pair : elevated_emis_var_names_) {
    const auto& var_name=pair.first;
    std::string item_name = "mam4_" + var_name + "_elevated_emiss_file_name";
    const auto file_name  = m_params.get<std::string>(item_name);
    util::TimeStamp ref_ts_vertical = mam_coupling::convert_date(elevated_emiss_cyclical_ymd);
    std::vector<Field> vertical_fields;
    for(const auto &field_name :pair.second) {
      vertical_fields.push_back(get_field_out(field_name+"_"+var_name).alias(field_name));
    }
    std::shared_ptr<DataInterpolation> di_vertical = std::make_shared<DataInterpolation>(grid_,vertical_fields);
    di_vertical->set_input_files_dimname(ShortFieldTagsNames::LEV,"altitude");
    di_vertical->setup_time_database ({file_name},util::TimeLine::YearlyPeriodic, ref_ts_vertical);
    di_vertical->create_horiz_remappers (extfrc_map_file=="none" ? "" : extfrc_map_file);
    DataInterpolation::VertRemapData remap_data_vertical;
    remap_data_vertical.vr_type = DataInterpolation::Custom;
    remap_data_vertical.pmid = z_iface;
    auto grid_after_hremap_vertical = di_vertical->get_grid_after_hremap();
    grid_after_hremap_vertical->reset_field_tag_name(ShortFieldTagsNames::LEV, "altitude");
    grid_after_hremap_vertical->reset_field_tag_name(ShortFieldTagsNames::ILEV, "altitude_int");
    // we create elevated emission remapper
    auto vertical_remapper_elevated = std::make_shared<VerticalRemapperMAM4>(grid_after_hremap_vertical, grid_,
    VerticalRemapperMAM4::VertRemapType::MAM4_ELEVATED_EMISSIONS);
    // we set source and target variables for interpolation
    vertical_remapper_elevated->set_source_pressure(file_name);
    vertical_remapper_elevated->set_target_pressure(z_iface);
    remap_data_vertical.custom_remapper=vertical_remapper_elevated;
    di_vertical->create_vert_remapper (remap_data_vertical);
    di_vertical->init_data_interval (start_of_step_ts());
    data_interp_elevated_emissions_.push_back(di_vertical);
  }//end var_name
}
#endif
// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMMicrophysics::initialize_impl(const RunType run_type) {
  // Determine orbital year. If orbital_year is negative, use current year
  // from timestamp for orbital year; if positive, use provided orbital year
  // for duration of simulation.
  // Check the interval values for the following fields used by this interface.
  // NOTE: We do not include aerosol and gas species, e.g., soa_a1, num_a1,
  // because we automatically added these fields.
  const std::map<std::string, std::pair<Real, Real>> ranges_microphysics = {
      {"SW_flux_dn", {-1e10, 1e10}},            // FIXME
      {"horiz_winds", {-1e10, 1e10}},           // FIXME
      {"nevapr", {-1e10, 1e10}},                // FIXME
      {"precip_ice_surf_mass", {-1e10, 1e10}},  // FIXME
      {"precip_liq_surf_mass", {-1e10, 1e10}},  // FIXME
      {"precip_total_tend", {-1e10, 1e10}},     // FIXME
      {"ps", {-1e10, 1e10}},                    // FIXME
      {"sfc_alb_dir_vis", {-1e10, 1e10}},       // FIXME
      {"snow_depth_land", {-1e10, 1e10}},       // FIXME
      {"surf_radiative_T", {-1e10, 1e10}},      // FIXME
      {"dqdt_h2so4_uptake", {-1e10, 1e10}}       // FIXME
  };
  set_ranges_process(ranges_microphysics);
  add_interval_checks();
  m_orbital_year = m_params.get<int>("orbital_year", -9999);

  // Get orbital parameters from yaml file
  m_orbital_eccen = m_params.get<double>("orbital_eccentricity", -9999);
  m_orbital_obliq = m_params.get<double>("orbital_obliquity", -9999);
  m_orbital_mvelp = m_params.get<double>("orbital_mvelp", -9999);

  // ---------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ---------------------------------------------------------------
  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);

  // FIXME: we are using cldfrac_tot in other mam4xx process.
  dry_atm_.cldfrac = get_field_in("cldfrac_liq").get_view<const Real **>();
  // FIXME: phis is not populated by populate_wet_and_dry_atm.
  dry_atm_.phis = get_field_in("phis").get_view<const Real *>();
  // get surface albedo: shortwave, direct
  d_sfc_alb_dir_vis_ = get_field_in("sfc_alb_dir_vis").get_view<const Real *>();

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

  if (extra_mam4_aero_microphys_diags_) {
    //Some dignostics fields have mixed units (kg/kg/s, #/kg/s, etc.)
    //For these fields, we add a docstring to the field to indicate that
    //the units are mixed and the user should be careful when using these fields.
    //Following map contains the map of fields with mixed units and their long names.
    const std::map<std::string, std::string>  mixed_units_fields = {
      {"mam4_microphysics_tendency_gas_phase_chemistry",
      "MAM4xx microphysics tendencies due to gas phase chemistry [mixed units: kg/kg/s or #/kg/s]"},

      {"mam4_microphysics_tendency_aqueous_chemistry",
      "MAM4xx microphysics tendencies due to aqueous chemistry [mixed units: kg/kg/s or #/kg/s]"},

      {"mam4_microphysics_tendency_aqso4",
      "MAM4xx microphysics tendencies due to aqueous SO4 [mixed units: kg/kg/s or #/kg/s]"},

      {"mam4_microphysics_tendency_aqh2so4",
      "MAM4xx microphysics tendencies due to aqueous H2SO4 [mixed units: kg/kg/s or #/kg/s]"},

      {"mam4_microphysics_tendency_condensation",
      "MAM4xx microphysics tendencies due to gas aerosol exchange (condensation) [mixed units: mol/mol/s or #/mol/s]"},

      {"mam4_microphysics_tendency_renaming",
      "MAM4xx microphysics tendencies due to gas aerosol exchange (renaming) [mixed units: mol/mol/s or #/mol/s]"},

      {"mam4_microphysics_tendency_nucleation",
      "MAM4xx microphysics tendencies due to gas aerosol exchange (nucleation) [mixed units: mol/mol/s or #/mol/s]"},

      {"mam4_microphysics_tendency_coagulation",
      "MAM4xx microphysics tendencies due to gas aerosol exchange (coagulation) [mixed units: mol/mol/s or #/mol/s]"},

      {"mam4_microphysics_tendency_renaming_cloud_borne",
      "MAM4xx microphysics tendencies due to gas aerosol exchange (renaming cloud borne) [mixed units: mol/mol/s or #/mol/s]"},

    };
    // Add docstring to the fields with mixed units
    add_io_docstring_to_fields_with_mixed_units(mixed_units_fields);
  }


  // set field property checks for the fields in this process
  /* e.g.
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,130.0,500.0,false);
  add_postcondition_check<LowerBound>(get_field_out("pbl_height"),m_grid,0);
  add_postcondition_check<Interval>(get_field_out("cldfrac_liq"),m_grid,0.0,1.0,false);
  add_postcondition_check<LowerBound>(get_field_out("tke"),m_grid,0);
  */

  if (config_.linoz.compute) {
    // climatology data for linear stratospheric chemistry
    auto linoz_o3_clim = buffer_.scratch[0];  // ozone (climatology) [vmr]
    auto linoz_o3col_clim =
        buffer_.scratch[1];  // column o3 above box (climatology) [Dobson Units
                             // (DU)]
    auto linoz_t_clim = buffer_.scratch[2];  // temperature (climatology) [K]
    auto linoz_PmL_clim =
        buffer_.scratch[3];  // P minus L (climatology) [vmr/s]
    auto linoz_dPmL_dO3 =
        buffer_.scratch[4];  // sensitivity of P minus L to O3 [1/s]
    auto linoz_dPmL_dT =
        buffer_.scratch[5];  // sensitivity of P minus L to T3 [K]
    auto linoz_dPmL_dO3col = buffer_.scratch[6];  // sensitivity of P minus L to
                                                  // overhead O3 column [vmr/DU]
    auto linoz_cariolle_pscs =
        buffer_.scratch[7];  // Cariolle parameter for PSC loss of ozone [1/s]
#endif
    auto ts = start_of_step_ts();
    std::string linoz_chlorine_file =
        m_params.get<std::string>("mam4_linoz_chlorine_file");
    int chlorine_loading_ymd = m_params.get<int>("mam4_chlorine_loading_ymd");
    scream::mam_coupling::create_linoz_chlorine_reader(
        linoz_chlorine_file, ts, chlorine_loading_ymd, chlorine_values_,
        chlorine_time_secs_);
  }

  init_temporary_views();
  // FIXME : why are we only using nlev_ instead of ncol_xnlev?
  cmfdqr_ = view_1d("cmfdqr_", nlev_);
  // Load the first month into extfrc_lst_end.
  // Note: At the first time step, the data will be moved into extfrc_lst_beg,
  //       and extfrc_lst_end will be reloaded from file with the new month.
  const int curr_month = start_of_step_ts().get_month() - 1;  // 0-based
#ifdef USE_OLD_LINOZ_FILE_READ
  if (config_.linoz.compute) {
    scream::mam_coupling::update_tracer_data_from_file(
      LinozDataReader_, curr_month, *LinozHorizInterp_, linoz_data_);
  }
  scream::mam_coupling::update_tracer_data_from_file(
      TracerDataReader_, curr_month, *TracerHorizInterp_, tracer_data_);
#else
  set_oxid_reader();
  set_linoz_reader();
#endif

#ifdef USE_OLD_VERTICAL_FILE_READ
  for(int i = 0; i < static_cast<int>(extfrc_lst_.size()); ++i) {
    scream::mam_coupling::update_tracer_data_from_file(
        ElevatedEmissionsDataReader_[i], curr_month,
        *ElevatedEmissionsHorizInterp_[i], elevated_emis_data_[i]);
  }
#else
  {
    set_elevated_emissions_reader();
    for(size_t i = 0; i < extfrc_lst_.size(); ++i) {
      std::string var_name = extfrc_lst_[i];
      const auto sector_names = elevated_emis_var_names_[var_name];
      const int nvars      = static_cast<int>(sector_names.size());
      forcings_[i].nsectors = nvars;
      // I am assuming the order of species in extfrc_lst_.
      // Indexing in mam4xx is fortran.
      forcings_[i].frc_ndx = i + 1;
      forcings_[i].file_alt_data = true;
      for(int isp = 0; isp < nvars; ++isp)
      {
        const std::string field_name = sector_names[isp]+"_"+var_name;
        forcings_[i].fields[isp] = get_field_out(field_name).get_view<Real **>();
      }//isp
    } //i
  }
#endif
  // //

  acos_cosine_zenith_host_ = view_1d_host("host_acos(cosine_zenith)", ncol_);
  acos_cosine_zenith_      = view_1d("device_acos(cosine_zenith)", ncol_);

}  // initialize_impl

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMMicrophysics::run_impl(const double dt) {
  const int ncol = ncol_;
  const int nlev = nlev_;
  //NOTE: get_default_team_policy produces a team size of 96 (nlev=72).
  // This interface hangs with this team size. Therefore,
  // let's use a team size of nlev.
#ifdef EKAT_ENABLE_GPU
       const int team_size=nlev;
#else
       const int team_size=1;
#endif
  const auto policy =
       ekat::ExeSpaceUtils<KT::ExeSpace>::get_team_policy_force_team_size(ncol, team_size);

  // preprocess input -- needs a scan for the calculation of atm height
  pre_process(wet_aero_, dry_aero_, wet_atm_, dry_atm_);
  Kokkos::fence();

  //----------- Variables from microphysics scheme -------------

  // Evaporation from stratiform rain [kg/kg/s]
  const auto &nevapr = get_field_in("nevapr").get_view<const Real **>();

  // Stratiform rain production rate [kg/kg/s]
  const auto &prain =
      get_field_in("precip_total_tend").get_view<const Real **>();

  const auto wet_geometric_mean_diameter_i =
      get_field_in("dgnumwet").get_view<const Real ***>();
  const auto dry_geometric_mean_diameter_i =
      get_field_in("dgnum").get_view<const Real ***>();
  const auto wetdens = get_field_in("wetdens").get_view<const Real ***>();

  // U wind component [m/s]
  const const_view_2d u_wind =
      get_field_in("horiz_winds").get_component(0).get_view<const Real **>();

  // V wind component [m/s]
  const const_view_2d v_wind =
      get_field_in("horiz_winds").get_component(1).get_view<const Real **>();

  // Liquid precip [kg/m2]
  const const_view_1d precip_liq_surf_mass =
      get_field_in("precip_liq_surf_mass").get_view<const Real *>();

  // Ice precip [kg/m2]
  const const_view_1d precip_ice_surf_mass =
      get_field_in("precip_ice_surf_mass").get_view<const Real *>();

  // Fractional land use [fraction]
  const const_view_2d fraction_landuse =
      get_field_in("fraction_landuse").get_view<const Real **>();

  // Downwelling solar flux at the surface [w/m2]
  const const_view_2d sw_flux_dn =
      get_field_in("SW_flux_dn").get_view<const Real **>();

  // Constituent fluxes of gas and aerosol species
  view_2d constituent_fluxes =
      get_field_out("constituent_fluxes").get_view<Real **>();

  // Surface temperature [K]
  const const_view_1d sfc_temperature =
      get_field_in("surf_radiative_T").get_view<const Real *>();

  // Surface pressure [Pa]
  const const_view_1d sfc_pressure =
      get_field_in("ps").get_view<const Real *>();

  // Snow depth on land [m]
  const const_view_1d snow_depth_land =
      get_field_in("snow_depth_land").get_view<const Real *>();

  // - dvmr/dt: Tendencies for mixing ratios  [kg/kg/s]
  view_2d dqdt_so4_aqueous_chemistry, dqdt_h2so4_uptake;
  view_3d gas_phase_chemistry_dvmrdt, aqueous_chemistry_dvmrdt;
  view_3d aqso4_incloud_mmr_tendency, aqh2so4_incloud_mmr_tendency;
  view_3d gas_aero_exchange_condensation, gas_aero_exchange_renaming,
          gas_aero_exchange_nucleation, gas_aero_exchange_coagulation,
          gas_aero_exchange_renaming_cloud_borne;

  if (extra_mam4_aero_microphys_diags_) {
    dqdt_so4_aqueous_chemistry = get_field_out("dqdt_so4_aqueous_chemistry").get_view<Real **>();
    dqdt_h2so4_uptake = get_field_out("dqdt_h2so4_uptake").get_view<Real **>();
    gas_phase_chemistry_dvmrdt = get_field_out("mam4_microphysics_tendency_gas_phase_chemistry").get_view<Real ***>();
    aqueous_chemistry_dvmrdt = get_field_out("mam4_microphysics_tendency_aqueous_chemistry").get_view<Real ***>();
    aqso4_incloud_mmr_tendency   = get_field_out("mam4_microphysics_tendency_aqso4").get_view<Real ***>();
    aqh2so4_incloud_mmr_tendency = get_field_out("mam4_microphysics_tendency_aqh2so4").get_view<Real ***>();
    gas_aero_exchange_condensation = get_field_out("mam4_microphysics_tendency_condensation").get_view<Real***>();
    gas_aero_exchange_renaming = get_field_out("mam4_microphysics_tendency_renaming").get_view<Real***>();
    gas_aero_exchange_nucleation = get_field_out("mam4_microphysics_tendency_nucleation").get_view<Real***>();
    gas_aero_exchange_coagulation = get_field_out("mam4_microphysics_tendency_coagulation").get_view<Real***>();
    gas_aero_exchange_renaming_cloud_borne = get_field_out("mam4_microphysics_tendency_renaming_cloud_borne").get_view<Real***>();
  }


  view_2d linoz_output[8];
  if (config_.linoz.compute) {
  // climatology data for linear stratospheric chemistry
  // ozone (climatology) [vmr]
  view_2d linoz_o3_clim;
  // column o3 above box (climatology) [Dobson Units (DU)]
  view_2d linoz_o3col_clim;
  // temperature (climatology) [K]
  view_2d linoz_t_clim;
  // P minus L (climatology) [vmr/s]
  view_2d linoz_PmL_clim;
  // sensitivity of P minus L to O3 [1/s]
  view_2d linoz_dPmL_dO3;
  // sensitivity of P minus L to T3 [K]
  view_2d linoz_dPmL_dT;
  // sensitivity of P minus L to overhead O3 column [vmr/DU]
  view_2d linoz_dPmL_dO3col;
  // Cariolle parameter for PSC loss of ozone [1/s]
  view_2d linoz_cariolle_pscs;


  linoz_output[0] = linoz_o3_clim;
  linoz_output[1] = linoz_o3col_clim;
  linoz_output[2] = linoz_t_clim;
  linoz_output[3] = linoz_PmL_clim;
  linoz_output[4] = linoz_dPmL_dO3;
  linoz_output[5] = linoz_dPmL_dT;
  linoz_output[6] = linoz_dPmL_dO3col;
  linoz_output[7] = linoz_cariolle_pscs;
  }
#else
  //FIXME: I will need to modify perform_atmospheric_chemistry_and_microphysics in MAM4xx.
  // Let's do that after we have tested the new interface for DataInterpolation.
  // const auto linoz_o3_clim = get_field_in("o3_clim").get_view<const Real **>();
  // const auto linoz_o3col_clim = get_field_in("o3col_clim").get_view<const Real **>();
  // const auto linoz_t_clim = get_field_in("t_clim").get_view<const Real **>();
  // const auto linoz_PmL_clim = get_field_in("PmL_clim").get_view<const Real **>();
  // const auto linoz_dPmL_dO3 = get_field_in("dPmL_dO3").get_view<const Real **>();
  // const auto linoz_dPmL_dT = get_field_in("dPmL_dT").get_view<const Real **>();
  // const auto linoz_dPmL_dO3col = get_field_in("dPmL_dO3col").get_view<const Real **>();
  // const auto linoz_cariolle_pscs = get_field_in("cariolle_pscs").get_view<const Real **>();
  data_interp_oxid_->run(end_of_step_ts());
  data_interp_linoz_->run(end_of_step_ts());

  const auto linoz_o3_clim = get_field_out("o3_clim").get_view<Real **>();
  const auto linoz_o3col_clim = get_field_out("o3col_clim").get_view<Real **>();
  const auto linoz_t_clim = get_field_out("t_clim").get_view<Real **>();
  const auto linoz_PmL_clim = get_field_out("PmL_clim").get_view<Real **>();
  const auto linoz_dPmL_dO3 = get_field_out("dPmL_dO3").get_view<Real **>();
  const auto linoz_dPmL_dT = get_field_out("dPmL_dT").get_view<Real **>();
  const auto linoz_dPmL_dO3col = get_field_out("dPmL_dO3col").get_view<Real **>();
  const auto linoz_cariolle_pscs = get_field_out("cariolle_pscs").get_view<Real **>();

  const auto oxid_O3 = get_field_out("O3").get_view<Real **>();
  const auto oxid_OH = get_field_out("OH").get_view<Real **>();
  const auto oxid_NO3 = get_field_out("NO3").get_view<Real **>();
  const auto oxid_HO2 = get_field_out("HO2").get_view<Real **>();

#endif
  // it's a bit wasteful to store this for all columns, but simpler from an
  // allocation perspective
  auto o3_col_dens = buffer_.scratch[8];

  /* Gather time and state information for interpolation */
  const auto ts = end_of_step_ts();

  if (config_.linoz.compute) {
  const Real chlorine_loading = scream::mam_coupling::chlorine_loading_advance(
      ts, chlorine_values_, chlorine_time_secs_);
    config_.linoz.chlorine_loading=chlorine_loading;
  }

#ifdef USE_OLD_LINOZ_FILE_READ
  // Update the TracerTimeState to reflect the current time
  trace_time_state_.t_now = ts.frac_of_year_in_days();
  scream::mam_coupling::advance_tracer_data(
      TracerDataReader_,                 // in
      *TracerHorizInterp_,               // out
      ts,                                // in
      trace_time_state_, tracer_data_,   // out
      dry_atm_.p_mid, dry_atm_.z_iface,  // in
      cnst_offline_);                    // out
  Kokkos::fence();
  const auto &cnst_offline = cnst_offline_;

  if (config_.linoz.compute) {
    view_2d linoz_output[8];
    linoz_output[0] = linoz_o3_clim;
    linoz_output[1] = linoz_o3col_clim;
    linoz_output[2] = linoz_t_clim;
    linoz_output[3] = linoz_PmL_clim;
    linoz_output[4] = linoz_dPmL_dO3;
    linoz_output[5] = linoz_dPmL_dT;
    linoz_output[6] = linoz_dPmL_dO3col;
    linoz_output[7] = linoz_cariolle_pscs;

    scream::mam_coupling::advance_tracer_data(
      LinozDataReader_,                  // in
      *LinozHorizInterp_,                // out
      ts,                                // in
      linoz_time_state_, linoz_data_,    // out
      dry_atm_.p_mid, dry_atm_.z_iface,  // in
      linoz_output);                     // out
    Kokkos::fence();
  }
#endif
#ifdef USE_OLD_VERTICAL_FILE_READ
  int i                            = 0;
  for(const auto &var_name : extfrc_lst_) {
    elevated_emiss_time_state_[i].t_now = ts.frac_of_year_in_days();
    const auto file_name = elevated_emis_file_name_[var_name];
    const auto var_names = elevated_emis_var_names_[var_name];
    auto& elevated_emis_output= forcings_[i].fields;
    scream::mam_coupling::advance_tracer_data(
        ElevatedEmissionsDataReader_[i], *ElevatedEmissionsHorizInterp_[i], ts,
        elevated_emiss_time_state_[i], elevated_emis_data_[i], dry_atm_.p_mid,
        dry_atm_.z_iface, elevated_emis_output);
    i++;
    Kokkos::fence();
  }
#else
  for (size_t i = 0; i < elevated_emis_var_names_.size(); ++i) {
    data_interp_elevated_emissions_[i]->run(end_of_step_ts());
  }
#endif

  const_view_1d &col_latitudes     = col_latitudes_;
  const_view_1d &d_sfc_alb_dir_vis = d_sfc_alb_dir_vis_;

  mam_coupling::DryAtmosphere &dry_atm = dry_atm_;
  mam_coupling::AerosolState &dry_aero = dry_aero_;


  mam4::mo_photo::PhotoTableData &photo_table = photo_table_;
  const Config &config                        = config_;
  const auto &work_photo_table                = work_photo_table_;
  const auto &photo_rates                     = photo_rates_;

  const auto &invariants   = invariants_;


  // Compute orbital parameters; these are used both for computing
  // the solar zenith angle.
  // Note: We are following the RRTMGP EAMxx interface to compute the zenith
  // angle. This operation is performed on the host because the routine
  // shr_orb_cosz_c2f has not been ported to C++.
  auto orbital_year = m_orbital_year;
  // Note: We need double precision because
  // shr_orb_params_c2f and shr_orb_decl_c2f only support double precision.
  double obliqr, lambm0, mvelpp;
  double eccen = m_orbital_eccen;
  double obliq = m_orbital_obliq;
  double mvelp = m_orbital_mvelp;
  // Use the orbital parameters to calculate the solar declination and
  // eccentricity factor
  double delta, eccf;
  if(eccen >= 0 && obliq >= 0 && mvelp >= 0) {
    // use fixed orbital parameters; to force this, we need to set
    // orbital_year to SHR_ORB_UNDEF_INT, which is exposed through
    // our c2f bridge as shr_orb_undef_int_c2f
    orbital_year = shr_orb_undef_int_c2f;
  } else if(orbital_year < 0) {
    // compute orbital parameters based on current year
    orbital_year = start_of_step_ts().get_year();
  }
  shr_orb_params_c2f(&orbital_year,                                       // in
                     &eccen, &obliq, &mvelp, &obliqr, &lambm0, &mvelpp);  // out

  // Want day + fraction; calday 1 == Jan 1 0Z
  auto calday = start_of_step_ts().frac_of_year_in_days() + 1;
  shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0, obliqr,  // in
                   &delta, &eccf);                         // out
  {
    const auto col_latitudes_host =
        grid_->get_geometry_data("lat").get_view<const Real *, Host>();
    const auto col_longitudes_host =
        grid_->get_geometry_data("lon").get_view<const Real *, Host>();
    // get a host copy of lat/lon
    // Determine the cosine zenith angle
    // NOTE: Since we are bridging to F90 arrays this must be done on HOST and
    // then deep copied to a device view.

    // Now use solar declination to calculate zenith angle for all points
    for(int i = 0; i < ncol; i++) {
      Real lat =
          col_latitudes_host(i) * M_PI / 180.0;  // Convert lat/lon to radians
      Real lon = col_longitudes_host(i) * M_PI / 180.0;
      // what's the aerosol microphys frequency?
      Real temp = shr_orb_cosz_c2f(calday, lat, lon, delta, dt);
      acos_cosine_zenith_host_(i) = acos(temp);
    }
    Kokkos::deep_copy(acos_cosine_zenith_, acos_cosine_zenith_host_);
  }
  const auto zenith_angle = acos_cosine_zenith_;
  constexpr int num_gas_aerosol_constituents = mam_coupling::gas_pcnst();

  const auto &extfrc   = extfrc_;
  const auto &forcings = forcings_;
  constexpr int extcnt = mam4::gas_chemistry::extcnt;

  const int offset_aerosol = mam4::utils::gasses_start_ind();
  Real adv_mass_kg_per_moles[num_gas_aerosol_constituents];
  // NOTE: Making copies of clsmap_4 and permute_4 to fix undefined arrays on
  // the device.
  int clsmap_4[num_gas_aerosol_constituents], permute_4[num_gas_aerosol_constituents];
  for(int i = 0; i < num_gas_aerosol_constituents; ++i) {
    // NOTE: state_q is kg/kg-dry-air; adv_mass is in g/mole.
    // Convert adv_mass to kg/mole as vmr_from_mmr function uses
    // molec_weight_dry_air with kg/mole units
    adv_mass_kg_per_moles[i] = mam4::gas_chemistry::adv_mass[i] / 1e3;
    clsmap_4[i]              = mam4::gas_chemistry::clsmap_4[i];
    permute_4[i]             = mam4::gas_chemistry::permute_4[i];
  }
  const auto &cmfdqr       = cmfdqr_;
  const auto &work_set_het = work_set_het_;
  const mam4::seq_drydep::Data drydep_data =
      mam4::seq_drydep::set_gas_drydep_data();
  const auto qv                = wet_atm_.qv;
  const int month              = start_of_step_ts().get_month();  // 1-based
  const int surface_lev        = nlev - 1;                 // Surface level
  const auto &index_season_lai = index_season_lai_;
  const int pcnst              = mam4::pcnst;
  const bool extra_mam4_aero_microphys_diags  = extra_mam4_aero_microphys_diags_;
  //NOTE: we need to initialize photo_rates_
  Kokkos::deep_copy(photo_rates_,0.0);
  // loop over atmosphere columns and compute aerosol microphysics

  Kokkos::parallel_for(
      "MAMMicrophysics::run_impl", policy,
      KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol     = team.league_rank();   // column index
        const Real col_lat = col_latitudes(icol);  // column latitude (degrees?)

        // convert column latitude to radians
        const Real rlats = col_lat * M_PI / 180.0;

        // fetch column-specific atmosphere state data
        const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
        const auto wet_diameter_icol =
            ekat::subview(wet_geometric_mean_diameter_i, icol);
        const auto dry_diameter_icol =
            ekat::subview(dry_geometric_mean_diameter_i, icol);
        const auto wetdens_icol = ekat::subview(wetdens, icol);

        // fetch column-specific subviews into aerosol prognostics
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);

        const auto invariants_icol = ekat::subview(invariants, icol);
        mam4::mo_setext::Forcing forcings_in[extcnt];

        for(int i = 0; i < extcnt; ++i) {
          const int nsectors       = forcings[i].nsectors;
          const int frc_ndx        = forcings[i].frc_ndx;
          const auto file_alt_data = forcings[i].file_alt_data;

          forcings_in[i].nsectors = nsectors;
          forcings_in[i].frc_ndx  = frc_ndx;
          // We may need to move this line where we read files.
          forcings_in[i].file_alt_data = file_alt_data;
          for(int isec = 0; isec < forcings[i].nsectors; ++isec) {
            const auto& field = forcings[i].fields[isec];
            forcings_in[i].fields_data[isec] = ekat::subview(field, icol);
          }
        }  // extcnt for loop

        const auto extfrc_icol = ekat::subview(extfrc, icol);

        view_1d cnst_offline_icol[mam4::mo_setinv::num_tracer_cnst];
 #ifdef USE_OLD_LINOZ_FILE_READ
        for(int i = 0; i < mam4::mo_setinv::num_tracer_cnst; ++i) {
          cnst_offline_icol[i] = ekat::subview(cnst_offline[i], icol);
        }
        // calculate o3 column densities (first component of col_dens in Fortran
        // code)
        auto o3_col_dens_i = ekat::subview(o3_col_dens, icol);
        const auto &work_photo_table_icol =
            ekat::subview(work_photo_table, icol);

        const auto &photo_rates_icol = ekat::subview(photo_rates, icol);

        mam4::microphysics::LinozData linoz_data;
        if (config.linoz.compute) {
          linoz_data.linoz_o3_clim_icol = ekat::subview(linoz_o3_clim, icol);
          linoz_data.linoz_t_clim_icol  = ekat::subview(linoz_t_clim, icol);
          linoz_data.linoz_o3col_clim_icol =
            ekat::subview(linoz_o3col_clim, icol);
          linoz_data.linoz_PmL_clim_icol = ekat::subview(linoz_PmL_clim, icol);
          linoz_data.linoz_dPmL_dO3_icol = ekat::subview(linoz_dPmL_dO3, icol);
          linoz_data.linoz_dPmL_dT_icol  = ekat::subview(linoz_dPmL_dT, icol);
          linoz_data.linoz_dPmL_dO3col_icol =
            ekat::subview(linoz_dPmL_dO3col, icol);
          linoz_data.linoz_cariolle_pscs_icol =
            ekat::subview(linoz_cariolle_pscs, icol);
        }
        const auto nevapr_icol       = ekat::subview(nevapr, icol);
        const auto prain_icol        = ekat::subview(prain, icol);
        const auto work_set_het_icol = ekat::subview(work_set_het, icol);

        mam4::MicrophysDiagnosticArrays diag_arrays;

        if (extra_mam4_aero_microphys_diags) {
          diag_arrays.dqdt_so4_aqueous_chemistry = ekat::subview(dqdt_so4_aqueous_chemistry, icol);
          diag_arrays.dqdt_h2so4_uptake          = ekat::subview(dqdt_h2so4_uptake, icol);
          diag_arrays.gas_phase_chemistry_dvmrdt = ekat::subview(gas_phase_chemistry_dvmrdt, icol);

          diag_arrays.aqueous_chemistry_dvmrdt   = ekat::subview(aqueous_chemistry_dvmrdt, icol);
          diag_arrays.aqso4_incloud_mmr_tendency = ekat::subview(aqso4_incloud_mmr_tendency, icol);
          diag_arrays.aqh2so4_incloud_mmr_tendency = ekat::subview(aqh2so4_incloud_mmr_tendency, icol);

          diag_arrays.gas_aero_exchange_condensation = ekat::subview(gas_aero_exchange_condensation, icol);
          diag_arrays.gas_aero_exchange_renaming = ekat::subview(gas_aero_exchange_renaming, icol);
          diag_arrays.gas_aero_exchange_nucleation = ekat::subview(gas_aero_exchange_nucleation, icol);
          diag_arrays.gas_aero_exchange_coagulation = ekat::subview(gas_aero_exchange_coagulation, icol);
          diag_arrays.gas_aero_exchange_renaming_cloud_borne = ekat::subview(gas_aero_exchange_renaming_cloud_borne, icol);
        }


        // Wind speed at the surface
        const Real wind_speed =
            haero::sqrt(u_wind(icol, surface_lev) * u_wind(icol, surface_lev) +
                        v_wind(icol, surface_lev) * v_wind(icol, surface_lev));

        // Total rain at the surface
        const Real rain =
            precip_liq_surf_mass(icol) + precip_ice_surf_mass(icol);

        // Snow depth on land [m]
        const Real snow_height = snow_depth_land(icol);

        // Downwelling solar flux at the surface (value at interface) [w/m2]
        const Real solar_flux = sw_flux_dn(icol, surface_lev + 1);

        Real fraction_landuse_icol[mam4::mo_drydep::n_land_type];
        for(int i = 0; i < mam4::mo_drydep::n_land_type; ++i) {
          fraction_landuse_icol[i] = fraction_landuse(icol, i);
        }
        int index_season[mam4::mo_drydep::n_land_type];
        {
          //-------------------------------------------------------------------------------------
          // define which season (relative to Northern hemisphere climate)
          //-------------------------------------------------------------------------------------

          //-------------------------------------------------------------------------------------
          // define season index based on fixed LAI
          //-------------------------------------------------------------------------------------
          for(int lt = 0; lt < mam4::mo_drydep::n_land_type; ++lt) {
            index_season[lt] = index_season_lai(icol, month - 1);
          }

          //-------------------------------------------------------------------------------------
          // special case for snow covered terrain
          //-------------------------------------------------------------------------------------
          if(snow_height > 0.01) {  // BAD_CONSTANT
            for(int lt = 0; lt < mam4::mo_drydep::n_land_type; ++lt) {
              index_season[lt] = 3;
            }
          }
        }
        // These output values need to be put somewhere:
        Real dflx_col[num_gas_aerosol_constituents] = {};  // deposition velocity [1/cm/s]
        Real dvel_col[num_gas_aerosol_constituents] = {};  // deposition flux [1/cm^2/s]
        // Output: values are dvel, dflx
        // Input/Output: progs::stateq, progs::qqcw
        team.team_barrier();
        mam4::microphysics::perform_atmospheric_chemistry_and_microphysics(
            team, dt, rlats, sfc_temperature(icol), sfc_pressure(icol),
            wind_speed, rain, solar_flux, cnst_offline_icol, forcings_in, atm,
            photo_table,  config.setsox, config.amicphys,
             zenith_angle(icol), d_sfc_alb_dir_vis(icol),
            o3_col_dens_i, photo_rates_icol, extfrc_icol, invariants_icol,
            work_photo_table_icol,
            config.linoz, linoz_data,
             eccf, adv_mass_kg_per_moles,
            fraction_landuse_icol, index_season, clsmap_4, permute_4,
            offset_aerosol,
            dry_diameter_icol, wet_diameter_icol,
            wetdens_icol, dry_atm.phis(icol), cmfdqr, prain_icol, nevapr_icol,
            work_set_het_icol, drydep_data, diag_arrays, dvel_col, dflx_col, progs);

        team.team_barrier();
        // Update constituent fluxes with gas drydep fluxes (dflx)
        // FIXME: Possible units mismatch (dflx is in kg/cm2/s but
        // constituent_fluxes is kg/m2/s) (Following mimics Fortran code
        // behavior but we should look into it)
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, offset_aerosol, pcnst), [&](int ispc) {
          constituent_fluxes(icol, ispc) -= dflx_col[ispc - offset_aerosol];
        });
      });  // parallel_for for the column loop
  Kokkos::fence();

  auto extfrc_fm = get_field_out("mam4_external_forcing").get_view<Real***>();

  // Avogadro's number [molecules/mol]
  const Real Avogadro = haero::Constants::avogadro;
  // Mapping from external forcing species index to physics constituent index
  // NOTE: These indices should match the species in extfrc_lst
  // TODO: getting rid of hard-coded indices
  Kokkos::Array<int, extcnt> extfrc_pcnst_index = {3, 6, 14, 27, 28, 13, 18, 30, 5};
  Kokkos::Array<Real, num_gas_aerosol_constituents> molar_mass_g_per_mol_tmp;
  for (int i = 0; i < num_gas_aerosol_constituents; ++i) {
    molar_mass_g_per_mol_tmp[i] = mam4::gas_chemistry::adv_mass[i];  // host-only access
  }

  // Transpose extfrc_ from internal layout [ncol][nlev][extcnt]
  // to output layout [ncol][extcnt][nlev]
  // This aligns with expected field storage in the EAMxx infrastructure.
  Kokkos::parallel_for("transpose_extfrc",
    Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0,0,0}, {ncol, extcnt, nlev}),
    KOKKOS_LAMBDA(const int i, const int j, const int k) {
      const int pcnst_idx = extfrc_pcnst_index[j];
      const Real molar_mass_g_per_mol = molar_mass_g_per_mol_tmp[pcnst_idx]; // g/mol
      // Modify units to MKS units: [molec/cm3/s] to [kg/m3/s]
      // Convert g → kg (× 1e-3), cm³ → m³ (× 1e6) → total factor: 1e-3 × 1e6 = 1e3 = 1000.0
      extfrc_fm(i,j,k) = extfrc(i,k,j) * (molar_mass_g_per_mol / Avogadro) * 1000.0;
  });

  // postprocess output
  post_process(wet_aero_, dry_aero_, dry_atm_);
  Kokkos::fence();

}  // MAMMicrophysics::run_impl

}  // namespace scream
