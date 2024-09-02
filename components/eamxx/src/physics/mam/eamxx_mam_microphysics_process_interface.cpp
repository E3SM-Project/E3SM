#include <netcdf.h>  // for serial NetCDF file reads on MPI root

#include <cmath>
#include <ekat/ekat_assert.hpp>
#include <physics/mam/eamxx_mam_microphysics_process_interface.hpp>
#include <share/io/scream_scorpio_interface.hpp>
#include <share/property_checks/field_lower_bound_check.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

#include "scream_config.h"  // for SCREAM_CIME_BUILD
#include "share/util/scream_common_physics_functions.hpp"

// NOTE: see the impl/ directory for the contents of the impl namespace
#include "impl/compute_o3_column_density.cpp"
#include "impl/gas_phase_chemistry.cpp"
#include "physics/rrtmgp/shr_orb_mod_c2f.hpp"

// For EKAT units package
#include "ekat/util/ekat_units.hpp"

// For EKAT subview
#include <ekat/kokkos/ekat_subview_utils.hpp>
#include <iostream>

namespace scream {

MAMMicrophysics::MAMMicrophysics(const ekat::Comm &comm,
                                 const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params), aero_config_() {
  configure(params);
}

AtmosphereProcessType MAMMicrophysics::type() const {
  return AtmosphereProcessType::Physics;
}

std::string MAMMicrophysics::name() const { return "mam4_micro"; }

namespace {

void set_data_file(const char *name, const char *path,
                   char location[MAX_FILENAME_LEN]) {
  EKAT_REQUIRE_MSG(strlen(SCREAM_DATA_DIR) + strlen(path) < MAX_FILENAME_LEN,
                   "Error! " << name << " path is too long (must be < "
                             << MAX_FILENAME_LEN << " characters)");
  sprintf(location, "%s/%s", SCREAM_DATA_DIR, path);
}

}  // namespace

#define set_file_location(data_file, path) \
  set_data_file(#data_file, path, data_file)

void MAMMicrophysics::set_defaults_() {
  config_.amicphys.do_cond   = true;
  config_.amicphys.do_rename = true;
  config_.amicphys.do_newnuc = true;
  config_.amicphys.do_coag   = true;

  config_.amicphys.nucleation                              = {};
  config_.amicphys.nucleation.dens_so4a_host               = 1770.0;
  config_.amicphys.nucleation.mw_so4a_host                 = 115.0;
  config_.amicphys.nucleation.newnuc_method_user_choice    = 2;
  config_.amicphys.nucleation.pbl_nuc_wang2008_user_choice = 1;
  config_.amicphys.nucleation.adjust_factor_pbl_ratenucl   = 1.0;
  config_.amicphys.nucleation.accom_coef_h2so4             = 1.0;
  config_.amicphys.nucleation.newnuc_adjust_factor_dnaitdt = 1.0;

  // these parameters guide the coupling between parameterizations
  // NOTE: mam4xx was ported with these parameters fixed, so it's probably not
  // NOTE: safe to change these without code modifications.
  config_.amicphys.gaexch_h2so4_uptake_optaa = 2;
  config_.amicphys.newnuc_h2so4_conc_optaa   = 2;

  //===========================================================
  // default data file locations (relative to SCREAM_DATA_DIR)
  //===========================================================

  // many of these paths were extracted from
  // e3smv2/bld/namelist_files/namelist_defaults_eam.xml

  // photolysis
  // set_file_location(config_.photolysis.rsf_file,
  // "../waccm/phot/RSF_GT200nm_v3.0_c080811.nc");
  // set_file_location(config_.photolysis.xs_long_file,
  // "../waccm/phot/temp_prs_GT200nm_JPL10_c130206.nc");

  // stratospheric chemistry
  // set_file_location(config_.linoz.chlorine_loading_file,
  // "../cam/chem/trop_mozart/ub/Linoz_Chlorine_Loading_CMIP6_0003-2017_c20171114.nc");

  // dry deposition
  set_file_location(config_.drydep.srf_file,
                    "../cam/chem/trop_mam/atmsrf_ne4pg2_200527.nc");
}

void MAMMicrophysics::configure(const ekat::ParameterList &params) {
  set_defaults_();
}

void MAMMicrophysics::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  Units nondim = Units::nondimensional();
  Units n_unit(1 / kg, "#/kg");  // number mixing ratios [# / kg air]
  const auto m2 = pow(m, 2);
  const auto s2 = pow(s, 2);

  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // number of levels per column

  // get column geometry and locations
  col_areas_      = grid_->get_geometry_data("area").get_view<const Real *>();
  col_latitudes_  = grid_->get_geometry_data("lat").get_view<const Real *>();
  col_longitudes_ = grid_->get_geometry_data("lon").get_view<const Real *>();

  // define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{{COL}, {ncol_}};

  // layout for 3D (2d horiz X 1d vertical) variables
  FieldLayout scalar3d_layout_mid{{COL, LEV}, {ncol_, nlev_}};
  // At interfaces
  FieldLayout scalar3d_layout_int{{COL, ILEV}, {ncol_, nlev_ + 1}};

  const int nmodes    = mam4::AeroConfig::num_modes();  // Number of modes
  // layout for 3D (ncol, nmodes, nlevs)
  FieldLayout scalar3d_mid_nmodes =
      grid_->get_3d_vector_layout(true, nmodes, "nmodes");

  static constexpr auto m3 = m * m * m;
     // Aerosol dry particle diameter [m]
  add_field<Required>("dgncur_a", scalar3d_mid_nmodes, m, grid_name);

  // Wet aerosol density [kg/m3]
  add_field<Required>("wetdens", scalar3d_mid_nmodes, kg / m3, grid_name);

  // Aerosol water [kg/kg]
  add_field<Required>("qaerwat", scalar3d_mid_nmodes, kg / kg, grid_name);

  // Wet Required diameter [m]
  add_field<Required>("dgnumwet", scalar3d_mid_nmodes, m, grid_name);

  // define fields needed in mam4xx

  // atmospheric quantities
  add_field<Required>("omega", scalar3d_layout_mid, Pa / s,
                      grid_name);  // vertical pressure velocity
  add_field<Required>("T_mid", scalar3d_layout_mid, K,
                      grid_name);  // Temperature
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa,
                      grid_name);  // total pressure
  // Total pressure [Pa] at interfaces
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid_name);
  add_field<Required>("qv", scalar3d_layout_mid, kg / kg, grid_name,
                      "tracers");  // specific humidity
  add_field<Required>("qi", scalar3d_layout_mid, kg / kg, grid_name,
                      "tracers");  // ice wet mixing ratio
  add_field<Required>("ni", scalar3d_layout_mid, n_unit, grid_name,
                      "tracers");  // ice number mixing ratio
  add_field<Required>("pbl_height", scalar2d_layout_col, m,
                      grid_name);  // planetary boundary layer height
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,
                      grid_name);  // p_del, hydrostatic pressure
  add_field<Required>("phis", scalar2d_layout_col, m2 / s2, grid_name);
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim,
                      grid_name);  // cloud fraction
  add_field<Required>("sfc_alb_dir_vis", scalar2d_layout_col, nondim,
                      grid_name);  // surface albedo shortwave, direct

  // droplet activation can alter cloud liquid and number mixing ratios
  add_field<Updated>("qc", scalar3d_layout_mid, kg / kg, grid_name,
                     "tracers");  // cloud liquid wet mixing ratio
  add_field<Updated>("nc", scalar3d_layout_mid, n_unit, grid_name,
                     "tracers");  // cloud liquid wet number mixing ratio

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name, "tracers");

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
    // NOT advected
    const char *cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);
    add_field<Updated>(cld_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name);

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);
      if(strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, kg / kg,
                           grid_name, "tracers");
      }

      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
      // NOT advected
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if(strlen(cld_mmr_field_name) > 0) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_layout_mid, kg / kg,
                           grid_name);
      }
    }  // end for loop for num species
  }    // end for loop for num modes

  // aerosol-related gases: mass mixing ratios
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_layout_mid, kg / kg,
                       grid_name, "tracers");
  }

  // Tracers group -- do we need this in addition to the tracers above? In any
  // case, this call should be idempotent, so it can't hurt.
  add_group<Updated>("tracers", grid_name, 1, Bundling::Required);
  // Creating a Linoz reader and setting Linoz parameters involves reading data
  // from a file and configuring the necessary parameters for the Linoz model.
  {
    linoz_file_name_ = m_params.get<std::string>("mam4_linoz_file_name");
    std::string spa_map_file = "";
    std::vector<std::string> var_names{"o3_clim",     "o3col_clim",   "t_clim",
                                       "PmL_clim",    "dPmL_dO3",     "dPmL_dT",
                                       "dPmL_dO3col", "cariolle_pscs"};

    // in format YYYYMMDD
    int linoz_cyclical_ymd = m_params.get<int>("mam4_linoz_ymd");
    scream::mam_coupling::setup_tracer_data(linoz_data_,
                           linoz_file_name_, linoz_cyclical_ymd);
    LinozHorizInterp_ = scream::mam_coupling::create_horiz_remapper(
        grid_, linoz_file_name_, spa_map_file, var_names, linoz_data_);
    LinozDataReader_ = scream::mam_coupling::create_tracer_data_reader(
        LinozHorizInterp_, linoz_file_name_);

    // linoz reader
    const auto io_grid_linoz = LinozHorizInterp_->get_src_grid();
    const int num_cols_io_linoz =
        io_grid_linoz->get_num_local_dofs();  // Number of columns on this rank
    const int num_levs_io_linoz =
        io_grid_linoz
            ->get_num_vertical_levels();  // Number of levels per column
    const int nvars = int(var_names.size());;
    linoz_data_.init(num_cols_io_linoz, num_levs_io_linoz, nvars);
    linoz_data_.allocate_temporal_views();
  }

  {
    oxid_file_name_          = m_params.get<std::string>("mam4_oxid_file_name");
    std::string spa_map_file = "";
    //NOTE: order matches mam4xx:
    std::vector<std::string> var_names{"O3", "OH", "NO3", "HO2"};

    // //in format YYYYMMDD
    int oxid_ymd = m_params.get<int>("mam4_oxid_ymd");
    scream::mam_coupling::setup_tracer_data(tracer_data_,
                           oxid_file_name_, oxid_ymd);
    TracerHorizInterp_ = scream::mam_coupling::create_horiz_remapper(
        grid_, oxid_file_name_, spa_map_file, var_names, tracer_data_);
    TracerDataReader_ = scream::mam_coupling::create_tracer_data_reader(
        TracerHorizInterp_, oxid_file_name_);

    const int nvars    = int(var_names.size());;
    const auto io_grid = TracerHorizInterp_->get_src_grid();
    const int num_cols_io =
        io_grid->get_num_local_dofs();  // Number of columns on this rank
    const int num_levs_io =
        io_grid->get_num_vertical_levels();  // Number of levels per column
    tracer_data_.init(num_cols_io, num_levs_io, nvars);
    tracer_data_.allocate_temporal_views();

    for(int ivar = 0; ivar < nvars; ++ivar) {
      cnst_offline_[ivar] = view_2d("cnst_offline_", ncol_, nlev_);
    }

  }

  {
    // FIXME: I will need to add this file per forcing file.
    std::string spa_map_file = "";
    // NOTE: order of forcing species is important.
    // extfrc_lst(:  9) = {'SO2             ','so4_a1          ','so4_a2          ','pom_a4          ','bc_a4           ',
                          // 'num_a1          ','num_a2          ','num_a4          ','SOAG            ' }
    // This order corresponds to files in namelist e3smv2
    extfrc_lst_=std::vector<std::string>({"so2","so4_a1","so4_a2","pom_a4","bc_a4",
                                         "num_a1","num_a2","num_a4","soag"});

    for (const auto& var_name : extfrc_lst_) {
      std::string item_name= "mam4_"+var_name+"_verti_emiss_file_name";
      const auto file_name = m_params.get<std::string>(item_name);
      vert_emis_file_name_[var_name] = file_name;
    }
    vert_emis_var_names_["so2"] = {"BB","ENE_ELEV", "IND_ELEV", "contvolc"};
    vert_emis_var_names_["so4_a1"] = {"BB","ENE_ELEV", "IND_ELEV", "contvolc"};
    vert_emis_var_names_["so4_a2"] = { "contvolc"};
    vert_emis_var_names_["pom_a4"] = {"BB"};
    vert_emis_var_names_["bc_a4"] = {"BB"};
    vert_emis_var_names_["num_a1"] = {"num_a1_SO4_ELEV_BB","num_a1_SO4_ELEV_ENE", "num_a1_SO4_ELEV_IND", "num_a1_SO4_ELEV_contvolc"};
    vert_emis_var_names_["num_a2"] = {"num_a2_SO4_ELEV_contvolc"};
    // num_a4
    // FIXME: why the sectors in this files are num_a1;
    //  I guess this should be num_a4? Is this a bug in the orginal nc files?
    vert_emis_var_names_["num_a4"] = {"num_a1_BC_ELEV_BB", "num_a1_POM_ELEV_BB"};
    vert_emis_var_names_["soag"] = {"SOAbb_src","SOAbg_src", "SOAff_src"};


    int verti_emiss_cyclical_ymd = m_params.get<int>("verti_emiss_ymd");

    for (const auto& var_name : extfrc_lst_) {
      const auto file_name = vert_emis_file_name_[var_name];
      const auto var_names = vert_emis_var_names_[var_name];

      scream::mam_coupling::TracerData data_tracer;
      scream::mam_coupling::setup_tracer_data(data_tracer,
                           file_name, verti_emiss_cyclical_ymd);
      auto hor_rem = scream::mam_coupling::create_horiz_remapper(
          grid_, file_name, spa_map_file, var_names, data_tracer);
      auto file_reader =
          scream::mam_coupling::create_tracer_data_reader(hor_rem, file_name);
      VertEmissionsHorizInterp_.push_back(hor_rem);
      VertEmissionsDataReader_.push_back(file_reader);
      vert_emis_data_.push_back(data_tracer);
    }// var_name vert emissions
    int i=0;
    int offset_emis_ver=0;
    for (auto it = extfrc_lst_.begin(); it != extfrc_lst_.end(); ++it, ++i) {
      const auto var_name = *it;
      const auto file_name = vert_emis_file_name_[var_name];
      const auto var_names = vert_emis_var_names_[var_name];
      const int nvars = int(var_names.size());

      forcings_[i].nsectors = nvars;
      // I am assuming the order of species in extfrc_lst_.
      // Indexing in mam4xx is fortran.
      forcings_[i].frc_ndx = i+1;
      const auto io_grid_emis = VertEmissionsHorizInterp_[i]->get_src_grid();
      const int num_cols_io_emis =
          io_grid_emis->get_num_local_dofs();  // Number of columns on this rank
      const int num_levs_io_emis =
          io_grid_emis
              ->get_num_vertical_levels();  // Number of levels per column
      vert_emis_data_[i].init(num_cols_io_emis, num_levs_io_emis, nvars);
      vert_emis_data_[i].allocate_temporal_views();
      forcings_[i].file_alt_data=vert_emis_data_[i].has_altitude_;
      for (int isp = 0; isp < nvars; ++isp)
      {
        EKAT_REQUIRE_MSG(
        offset_emis_ver <= int(mam_coupling::MAX_NUM_VERT_EMISSION_FIELDS),
        "Error! Number of fields is bigger than MAX_NUM_VERT_EMISSION_FIELDS. Increase the MAX_NUM_VERT_EMISSION_FIELDS in helper_micro.hpp \n");
        forcings_[i].offset=offset_emis_ver;
        vert_emis_output_[isp+offset_emis_ver] =
          view_2d("vert_emis_output_", ncol_, nlev_);
      }
      offset_emis_ver+=nvars;
    }  // end i

  }
}

// this checks whether we have the tracers we expect
void MAMMicrophysics::set_computed_group_impl(const FieldGroup &group) {
  const auto &name = group.m_info->m_group_name;
  EKAT_REQUIRE_MSG(
      name == "tracers",
      "Error! MAM4 expects a 'tracers' field group (got '" << name << "')\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
                   "Error! MAM4 expects bundled fields for tracers.\n");

  // how many aerosol/gas tracers do we expect?
  int num_tracers = mam_coupling::num_aero_modes() +
                    mam_coupling::num_aero_tracers() +
                    mam_coupling::num_aero_gases();
  EKAT_REQUIRE_MSG(group.m_info->size() >= num_tracers,
                   "Error! MAM4 requires at least "
                       << group.m_info->size() << " " << num_tracers << " "
                       << mam_coupling::num_aero_modes() << " "
                       << mam_coupling::num_aero_tracers() << " "
                       << mam_coupling::num_aero_gases() << "aerosol tracers.");
}

size_t MAMMicrophysics::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

void MAMMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(
      used_mem == requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for MAMMicrophysics.");
}

void MAMMicrophysics::initialize_impl(const RunType run_type) {
  step_ = 0;

  // Determine orbital year. If orbital_year is negative, use current year
  // from timestamp for orbital year; if positive, use provided orbital year
  // for duration of simulation.
  m_orbital_year = m_params.get<Int>("orbital_year", -9999);

  // Get orbital parameters from yaml file
  m_orbital_eccen = m_params.get<Int>("orbital_eccentricity", -9999);
  m_orbital_obliq = m_params.get<Int>("orbital_obliquity", -9999);
  m_orbital_mvelp = m_params.get<Int>("orbital_mvelp", -9999);

  // populate the wet and dry atmosphere states with views from fields and
  // the buffer
  wet_atm_.qv = get_field_in("qv").get_view<const Real **>();
  wet_atm_.qc = get_field_out("qc").get_view<Real **>();
  wet_atm_.nc = get_field_out("nc").get_view<Real **>();
  wet_atm_.qi = get_field_in("qi").get_view<const Real **>();
  wet_atm_.ni = get_field_in("ni").get_view<const Real **>();

  dry_atm_.T_mid   = get_field_in("T_mid").get_view<const Real **>();
  dry_atm_.p_mid   = get_field_in("p_mid").get_view<const Real **>();
  dry_atm_.p_int   = get_field_in("p_int").get_view<const Real **>();
  dry_atm_.p_del   = get_field_in("pseudo_density").get_view<const Real **>();
  dry_atm_.cldfrac = get_field_in("cldfrac_tot")
                         .get_view<const Real **>();  // FIXME: tot or liq?
  dry_atm_.pblh      = get_field_in("pbl_height").get_view<const Real *>();
  dry_atm_.phis      = get_field_in("phis").get_view<const Real *>();
  dry_atm_.omega     = get_field_in("omega").get_view<const Real **>();
  dry_atm_.z_mid     = buffer_.z_mid;
  dry_atm_.dz        = buffer_.dz;
  dry_atm_.z_iface   = buffer_.z_iface;
  dry_atm_.qv        = buffer_.qv_dry;
  dry_atm_.qc        = buffer_.qc_dry;
  dry_atm_.nc        = buffer_.nc_dry;
  dry_atm_.qi        = buffer_.qi_dry;
  dry_atm_.ni        = buffer_.ni_dry;
  dry_atm_.w_updraft = buffer_.w_updraft;
  dry_atm_.z_surf    = 0.0;  // It is always zero.

  // get surface albedo: shortwave, direct
  d_sfc_alb_dir_vis_ = get_field_in("sfc_alb_dir_vis").get_view<const Real *>();

  // perform any initialization work
  if(run_type == RunType::Initial) {
  }

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    wet_aero_.int_aero_nmr[m] =
        get_field_out(int_nmr_field_name).get_view<Real **>();
    dry_aero_.int_aero_nmr[m] = buffer_.dry_int_aero_nmr[m];

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    const char *cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);
    wet_aero_.cld_aero_nmr[m] =
        get_field_out(cld_nmr_field_name).get_view<Real **>();
    dry_aero_.cld_aero_nmr[m] = buffer_.dry_cld_aero_nmr[m];

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);
      if(strlen(int_mmr_field_name) > 0) {
        wet_aero_.int_aero_mmr[m][a] =
            get_field_out(int_mmr_field_name).get_view<Real **>();
        dry_aero_.int_aero_mmr[m][a] = buffer_.dry_int_aero_mmr[m][a];
      }

      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if(strlen(cld_mmr_field_name) > 0) {
        wet_aero_.cld_aero_mmr[m][a] =
            get_field_out(cld_mmr_field_name).get_view<Real **>();
        dry_aero_.cld_aero_mmr[m][a] = buffer_.dry_cld_aero_mmr[m][a];
      }
    }
  }

  // set wet/dry aerosol-related gas state data
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] = get_field_out(mmr_field_name).get_view<Real **>();
    dry_aero_.gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

  // create our photolysis rate calculation table
  const std::string rsf_file = m_params.get<std::string>("mam4_rsf_file");
  const std::string xs_long_file =
      m_params.get<std::string>("mam4_xs_long_file");

  photo_table_ = impl::read_photo_table(rsf_file, xs_long_file);
  // set up our preprocess/postprocess functors
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                         dry_aero_);
  postprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                          dry_aero_);

  // set field property checks for the fields in this process
  /* e.g.
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,130.0,500.0,false);
  add_postcondition_check<LowerBound>(get_field_out("pbl_height"),m_grid,0);
  add_postcondition_check<Interval>(get_field_out("cldfrac_liq"),m_grid,0.0,1.0,false);
  add_postcondition_check<LowerBound>(get_field_out("tke"),m_grid,0);
  */

  // set up WSM for internal local variables
  // (we'll probably need this later, but we'll just use ATMBufferManager for
  // now)
  // const auto default_policy =
  // ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  // workspace_mgr_.setup(buffer_.wsm_data, nlev_+1,
  // 13+(n_wind_slots+n_trac_slots), default_policy);

  {
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

    // const std::string linoz_chlorine_file =
    // "Linoz_Chlorine_Loading_CMIP6_0003-2017_c20171114.nc"; auto ts =
    // timestamp(); int chlorine_loading_ymd=20100101;
    auto ts = timestamp();
    std::string linoz_chlorine_file =
        m_params.get<std::string>("mam4_linoz_chlorine_file");
    int chlorine_loading_ymd = m_params.get<int>("mam4_chlorine_loading_ymd");
    scream::mam_coupling::create_linoz_chlorine_reader(
        linoz_chlorine_file, ts, chlorine_loading_ymd, chlorine_values_,
        chlorine_time_secs_);
  }

  const int photo_table_len = get_photo_table_work_len(photo_table_);
  work_photo_table_ = view_2d("work_photo_table", ncol_, photo_table_len);

  // here's where we store per-column photolysis rates
  photo_rates_ = view_3d("photo_rates", ncol_, nlev_, mam4::mo_photo::phtcnt);

  //
  // Load the first month into spa_end.
  // Note: At the first time step, the data will be moved into spa_beg,
  //       and spa_end will be reloaded from file with the new month.
  const int curr_month = timestamp().get_month() - 1;  // 0-based

  scream::mam_coupling::update_tracer_data_from_file(
        LinozDataReader_,
        curr_month,
        *LinozHorizInterp_,
        linoz_data_);

  scream::mam_coupling::update_tracer_data_from_file(
        TracerDataReader_, curr_month, *TracerHorizInterp_,
        tracer_data_);

  int i=0;
  for (auto it = extfrc_lst_.begin(); it != extfrc_lst_.end(); ++it, ++i) {
    scream::mam_coupling::update_tracer_data_from_file(
          VertEmissionsDataReader_[i], curr_month,
          *VertEmissionsHorizInterp_[i], vert_emis_data_[i]);
  }

  invariants_ = view_3d("invarians", ncol_, nlev_, mam4::gas_chemistry::nfs);

  constexpr int extcnt =
      mam4::gas_chemistry::extcnt;
  extfrc_=view_3d("extfrc_", ncol_, nlev_, extcnt);
}

void MAMMicrophysics::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);
  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  const auto wet_geometric_mean_diameter_i = get_field_in("dgnumwet").get_view<const Real ***>();
  const auto dry_geometric_mean_diameter_i = get_field_in("dgncur_a").get_view<const  Real ***>();
  const auto qaerwat = get_field_in("qaerwat").get_view<const Real ***>();
  const auto wetdens = get_field_in("wetdens").get_view<const Real ***>();

  // reset internal WSM variables
  // workspace_mgr_.reset_internals();

  // NOTE: nothing depends on simulation time (yet), so we can just use zero for
  // now
  // double t = 0.0;

  // climatology data for linear stratospheric chemistry
  auto linoz_o3_clim = buffer_.scratch[0];  // ozone (climatology) [vmr]
  auto linoz_o3col_clim =
      buffer_
          .scratch[1];  // column o3 above box (climatology) [Dobson Units (DU)]
  auto linoz_t_clim   = buffer_.scratch[2];  // temperature (climatology) [K]
  auto linoz_PmL_clim = buffer_.scratch[3];  // P minus L (climatology) [vmr/s]
  auto linoz_dPmL_dO3 =
      buffer_.scratch[4];  // sensitivity of P minus L to O3 [1/s]
  auto linoz_dPmL_dT =
      buffer_.scratch[5];  // sensitivity of P minus L to T3 [K]
  auto linoz_dPmL_dO3col = buffer_.scratch[6];  // sensitivity of P minus L to
                                                // overhead O3 column [vmr/DU]
  auto linoz_cariolle_pscs =
      buffer_.scratch[7];  // Cariolle parameter for PSC loss of ozone [1/s]

  view_2d linoz_output[8];
  linoz_output[0] = linoz_o3_clim;
  linoz_output[1] = linoz_o3col_clim;
  linoz_output[2] = linoz_t_clim;
  linoz_output[3] = linoz_PmL_clim;
  linoz_output[4] = linoz_dPmL_dO3;
  linoz_output[5] = linoz_dPmL_dT;
  linoz_output[6] = linoz_dPmL_dO3col;
  linoz_output[7] = linoz_cariolle_pscs;

  // it's a bit wasteful to store this for all columns, but simpler from an
  // allocation perspective
  auto o3_col_dens = buffer_.scratch[8];

  /* Gather time and state information for interpolation */
  const auto ts = timestamp() + dt;

  const Real chlorine_loading = scream::mam_coupling::chlorine_loading_advance(
      ts, chlorine_values_, chlorine_time_secs_);

  // /* Update the TracerTimeState to reflect the current time, note the
  // addition of dt */
  trace_time_state_.t_now = ts.frac_of_year_in_days();
  scream::mam_coupling::advance_tracer_data(
      TracerDataReader_, *TracerHorizInterp_, ts, trace_time_state_,
      tracer_data_,
      dry_atm_.p_mid, dry_atm_.z_iface, cnst_offline_);

  scream::mam_coupling::advance_tracer_data(
      LinozDataReader_, *LinozHorizInterp_, ts, linoz_time_state_,
      linoz_data_,
      dry_atm_.p_mid, dry_atm_.z_iface, linoz_output);

  vert_emiss_time_state_.t_now = ts.frac_of_year_in_days();
  int i=0;
  for (auto it = extfrc_lst_.begin(); it != extfrc_lst_.end(); ++it, ++i) {
    const auto var_name = *it;
    const auto file_name = vert_emis_file_name_[var_name];
    const auto var_names = vert_emis_var_names_[var_name];
    const int nsectors = int(var_names.size());
    view_2d vert_emis_output[nsectors];
    for (int isp = 0; isp < nsectors; ++isp)
    {
      vert_emis_output[isp]= vert_emis_output_[isp+forcings_[i].offset];
    }
    scream::mam_coupling::advance_tracer_data(
        VertEmissionsDataReader_[i], *VertEmissionsHorizInterp_[i], ts,
        vert_emiss_time_state_,
        vert_emis_data_[i], dry_atm_.p_mid,
        dry_atm_.z_iface, vert_emis_output);
  }

  const_view_1d &col_latitudes     = col_latitudes_;
  const_view_1d &col_longitudes    = col_longitudes_;
  const_view_1d &d_sfc_alb_dir_vis = d_sfc_alb_dir_vis_;



  mam_coupling::DryAtmosphere &dry_atm        = dry_atm_;
  mam_coupling::AerosolState &dry_aero        = dry_aero_;
  mam4::mo_photo::PhotoTableData &photo_table = photo_table_;
  const int nlev                              = nlev_;
  const Config &config                        = config_;
  const auto &step                            = step_;
  const auto &work_photo_table = work_photo_table_;
  const auto &photo_rates      = photo_rates_;

  const auto &invariants   = invariants_;
  const auto &cnst_offline = cnst_offline_;

  // Compute orbital parameters; these are used both for computing
  // the solar zenith angle.
  auto ts2 = timestamp();
  double obliqr, lambm0, mvelpp;
  auto orbital_year = m_orbital_year;
  auto eccen        = m_orbital_eccen;
  auto obliq        = m_orbital_obliq;
  auto mvelp        = m_orbital_mvelp;
  if(eccen >= 0 && obliq >= 0 && mvelp >= 0) {
    // use fixed orbital parameters; to force this, we need to set
    // orbital_year to SHR_ORB_UNDEF_INT, which is exposed through
    // our c2f bridge as shr_orb_undef_int_c2f
    orbital_year = shr_orb_undef_int_c2f;
  } else if(orbital_year < 0) {
    // compute orbital parameters based on current year
    orbital_year = ts2.get_year();
  }
  shr_orb_params_c2f(&orbital_year, &eccen, &obliq, &mvelp, &obliqr, &lambm0,
                     &mvelpp);
  // Use the orbital parameters to calculate the solar declination and
  // eccentricity factor
  Real delta, eccf;
  auto calday = ts2.frac_of_year_in_days() +
                1;  // Want day + fraction; calday 1 == Jan 1 0Z
  shr_orb_decl_c2f(calday, eccen, mvelpp, lambm0, obliqr, &delta, &eccf);
  constexpr int num_modes = mam4::AeroConfig::num_modes();
  constexpr int gas_pcnst = mam_coupling::gas_pcnst();
  constexpr int nqtendbb  = mam_coupling::nqtendbb();

  constexpr int pcnst = mam4::pcnst;
  const auto vert_emis_output = vert_emis_output_;
  const auto extfrc = extfrc_;
  const auto forcings = forcings_;
  constexpr int extcnt = mam4::gas_chemistry::extcnt;

  // FIXME: remove this hard-code value
  const int offset_aerosol = mam4::utils::gasses_start_ind();
  Real adv_mass_kg_per_moles[gas_pcnst];
  // NOTE: Making copies of clsmap_4 and permute_4 to fix undefined arrays on the device.
  int clsmap_4[gas_pcnst], permute_4[gas_pcnst];
  for(int i = 0; i < gas_pcnst; ++i)
  {
    adv_mass_kg_per_moles[i] = mam4::gas_chemistry::adv_mass[i]/1e3;
    clsmap_4[i]=mam4::gas_chemistry::clsmap_4[i];
    permute_4[i]=mam4::gas_chemistry::permute_4[i];
  }

  // loop over atmosphere columns and compute aerosol microphyscs
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol = team.league_rank();  // column index

        Real col_lat = col_latitudes(icol);   // column latitude (degrees?)
        Real col_lon = col_longitudes(icol);  // column longitude

        Real rlats =
            col_lat * M_PI / 180.0;  // convert column latitude to radians
        Real rlons =
            col_lon * M_PI / 180.0;  // convert column longitude to radians

        // fetch column-specific atmosphere state data
        auto atm     = mam_coupling::atmosphere_for_column(dry_atm, icol);
        auto z_iface = ekat::subview(dry_atm.z_iface, icol);
        Real phis    = dry_atm.phis(icol);

        auto wet_diameter_icol =
            ekat::subview(wet_geometric_mean_diameter_i, icol);
        auto dry_diameter_icol =
            ekat::subview(dry_geometric_mean_diameter_i, icol);
        auto qaerwat_icol     = ekat::subview(qaerwat, icol);
        auto wetdens_icol     = ekat::subview(wetdens, icol);


        // set surface state data
        haero::Surface sfc{};

        // fetch column-specific subviews into aerosol prognostics
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);

        // set up diagnostics
        mam4::Diagnostics diags(nlev);
        const auto invariants_icol = ekat::subview(invariants, icol);
        mam4::mo_setext::Forcing forcings_in[extcnt];
        for (int i = 0; i < extcnt; ++i)
        {
          int nsectors = forcings[i].nsectors;
          int frc_ndx = forcings[i].frc_ndx;
          auto file_alt_data= forcings[i].file_alt_data;

          forcings_in[i].nsectors=nsectors;
	        forcings_in[i].frc_ndx=frc_ndx;
          // We may need to move this line where we read files.
          forcings_in[i].file_alt_data=file_alt_data;
          for (int isec = 0; isec < forcings[i].nsectors; ++isec)
          {
          const auto field = vert_emis_output[isec+forcings[i].offset];
          forcings_in[i].fields_data[isec]=ekat::subview(field, icol);
          }
        }
        const auto extfrc_icol = ekat::subview(extfrc, icol);

        mam4::mo_setext::extfrc_set(forcings_in, extfrc_icol);

        view_1d cnst_offline_icol[mam4::mo_setinv::num_tracer_cnst];
        for (int i = 0; i < mam4::mo_setinv::num_tracer_cnst; ++i) {
            cnst_offline_icol[i] = ekat::subview(cnst_offline[i],icol);
        }

        mam4::mo_setinv::setinv(team, invariants_icol, atm.temperature,
                                atm.vapor_mixing_ratio,
                                cnst_offline_icol, atm.pressure);

        // calculate o3 column densities (first component of col_dens in Fortran
        // code)
        auto o3_col_dens_i = ekat::subview(o3_col_dens, icol);
        impl::compute_o3_column_density(team, atm, progs, invariants_icol,
                                        o3_col_dens_i);

        // set up photolysis work arrays for this column.
        mam4::mo_photo::PhotoTableWorkArrays photo_work_arrays_icol;
        const auto& work_photo_table_icol = ekat::subview(work_photo_table,
        icol);
        //  set work view using 1D photo_work_arrays_icol
        mam4::mo_photo::set_photo_table_work_arrays(photo_table,
                                                    work_photo_table_icol,
                                                    photo_work_arrays_icol);

        // ... look up photolysis rates from our table
        // NOTE: the table interpolation operates on an entire column of data,
        // so we NOTE: must do it before dispatching to individual vertical
        // levels
        Real zenith_angle =
            shr_orb_cosz_c2f(calday, rlats, rlons, delta,
                             dt);  // what's the aerosol microphys frequency?
        zenith_angle = acos(zenith_angle);

        Real surf_albedo = d_sfc_alb_dir_vis(icol);

        const auto &photo_rates_icol = ekat::subview(photo_rates, icol);

        mam4::mo_photo::table_photo(photo_rates_icol, atm.pressure,
        atm.hydrostatic_dp,
         atm.temperature, o3_col_dens_i, zenith_angle, surf_albedo,
         atm.liquid_mixing_ratio, atm.cloud_fraction, eccf, photo_table,
         photo_work_arrays_icol);

        // compute aerosol microphysics on each vertical level within this
        // column
        Kokkos::parallel_for(
            Kokkos::TeamThreadRange(team, nlev), [&](const int k) {
              // extract atm state variables (input)
              Real temp    = atm.temperature(k);
              Real pmid    = atm.pressure(k);
              Real pdel    = atm.hydrostatic_dp(k);
              Real zm      = atm.height(k);
              Real zi      = z_iface(k);
              Real pblh    = atm.planetary_boundary_layer_height;
              Real qv      = atm.vapor_mixing_ratio(k);
              Real cldfrac = atm.cloud_fraction(k);
              Real lwc     = atm.liquid_mixing_ratio(k);
              Real cldnum  = atm.cloud_liquid_number_mixing_ratio(k);

              // extract aerosol state variables into "working arrays" (mass
              // mixing ratios) (in EAM, this is done in the gas_phase_chemdr
              // subroutine defined within
              //  mozart/mo_gas_phase_chemdr.F90)
              Real q[gas_pcnst]    = {};
              Real qqcw[gas_pcnst] = {};
              // mam_coupling::transfer_prognostics_to_work_arrays(progs, k, q,
              //                                                   qqcw);
              Real state_q[pcnst]    = {};
              Real qqcw_long[pcnst] = {};
              mam4::utils::extract_stateq_from_prognostics(progs,atm, state_q, k);

              mam4::utils::extract_qqcw_from_prognostics(progs,qqcw_long,k);

              for (int i = offset_aerosol; i < pcnst; ++i) {
                q[i-offset_aerosol] =state_q[i];
                qqcw[i-offset_aerosol] =qqcw_long[i];
              }
              // convert mass mixing ratios to volume mixing ratios (VMR),
              // equivalent to tracer mixing ratios (TMR))
              Real vmr[gas_pcnst], vmrcw[gas_pcnst];
              // CHECK: convert_work_arrays_to_vmr and mmr2vmr should produce the same ouputs
              // However, in mmr2vmr we do not iterate to get species value.
              // mam_coupling::convert_work_arrays_to_vmr(q, qqcw, vmr, vmrcw);
              mam_coupling::mmr2vmr(q,adv_mass_kg_per_moles, vmr);
              mam_coupling::mmr2vmr(qqcw,adv_mass_kg_per_moles,vmrcw);

              // aerosol/gas species tendencies (output)
              Real vmr_tendbb[gas_pcnst][nqtendbb]   = {};
              Real vmrcw_tendbb[gas_pcnst][nqtendbb] = {};

              // create work array copies to retain "pre-chemistry" values
              Real vmr_pregaschem[gas_pcnst]   = {};
              Real vmr_precldchem[gas_pcnst]   = {};
              Real vmrcw_precldchem[gas_pcnst] = {};
              for(int i = 0; i < gas_pcnst; ++i) {
                vmr_pregaschem[i]   = vmr[i];
                vmr_precldchem[i]   = vmr[i];
                vmrcw_precldchem[i] = vmrcw[i];
              }

              //---------------------
              // Gas Phase Chemistry
              //---------------------
              //
              const auto& extfrc_k = ekat::subview(extfrc_icol,k);
              const auto& invariants_k = ekat::subview(invariants_icol,k);
              const auto& photo_rates_k = ekat::subview(photo_rates_icol,k);
              impl::gas_phase_chemistry(zm, zi, phis, temp, pmid, pdel, dt,
                                        photo_rates_k.data(), extfrc_k.data(),
                                        invariants_k.data(), clsmap_4, permute_4, vmr);
              //----------------------
              // Aerosol microphysics
              //----------------------
              // the logic below is taken from the aero_model_gasaerexch
              // subroutine in eam/src/chemistry/modal_aero/aero_model.F90

              // aqueous chemistry ...
              const int loffset =
                  8;  // offset of first tracer in work arrays
                      // (taken from mam4xx setsox validation test)
              const Real mbar      = haero::Constants::molec_weight_dry_air;
              constexpr int indexm = mam4::gas_chemistry::indexm;
              mam4::mo_setsox::setsox_single_level(
                  loffset, dt, pmid, pdel, temp, mbar, lwc, cldfrac, cldnum,
                  invariants_k[indexm], config.setsox, vmrcw, vmr);

              // calculate aerosol water content using water uptake treatment
              // * dry and wet diameters [m]
              // * wet densities [kg/m3]
              // * aerosol water mass mixing ratio [kg/kg]
              Real dgncur_a_kk[num_modes]    = {};
              Real dgncur_awet_kk[num_modes] = {};
              Real wetdens_kk[num_modes]     = {};
              Real qaerwat_kk[num_modes]     = {};

              for (int imode = 0; imode < num_modes; imode++) {
                dgncur_awet_kk[imode] = wet_diameter_icol(imode, k);
                dgncur_a_kk[imode] = dry_diameter_icol(imode, k);
                qaerwat_kk[imode] = qaerwat_icol(imode, k);
                wetdens_kk[imode] = wetdens_icol(imode, k);
              }
              // do aerosol microphysics (gas-aerosol exchange, nucleation,
              // coagulation)
              impl::modal_aero_amicphys_intr(
                  config.amicphys, step, dt, temp, pmid, pdel, zm, pblh, qv,
                  cldfrac, vmr, vmrcw, vmr_pregaschem, vmr_precldchem,
                  vmrcw_precldchem, vmr_tendbb, vmrcw_tendbb, dgncur_a_kk,
                  dgncur_awet_kk, wetdens_kk, qaerwat_kk);
              //-----------------
              // LINOZ chemistry
              //-----------------

              // the following things are diagnostics, which we're not
              // including in the first rev
              Real do3_linoz, do3_linoz_psc, ss_o3, o3col_du_diag,
                  o3clim_linoz_diag, zenith_angle_degrees;

              int o3_ndx = 0;  // index of "O3" in solsym array (in EAM)

      mam4::lin_strat_chem::lin_strat_chem_solve_kk(o3_col_dens_i(k), temp,
        zenith_angle, pmid, dt, rlats,
        linoz_o3_clim(icol, k), linoz_t_clim(icol, k), linoz_o3col_clim(icol, k),
        linoz_PmL_clim(icol, k), linoz_dPmL_dO3(icol, k), linoz_dPmL_dT(icol, k),
        linoz_dPmL_dO3col(icol, k), linoz_cariolle_pscs(icol, k),
        chlorine_loading, config.linoz.psc_T, vmr[o3_ndx],
        do3_linoz, do3_linoz_psc, ss_o3,
        o3col_du_diag, o3clim_linoz_diag, zenith_angle_degrees);

              // update source terms above the ozone decay threshold
              /*if (k > nlev - config.linoz.o3_lbl - 1) {
                Real do3_mass; // diagnostic, not needed
                mam4::lin_strat_chem::lin_strat_sfcsink_kk(dt, pdel,
              vmr[o3_ndx], config.linoz.o3_sfc, config.linoz.o3_tau, do3_mass);
              }*/

              // ... check for negative values and reset to zero
              for(int i = 0; i < gas_pcnst; ++i) {
                if(vmr[i] < 0.0) vmr[i] = 0.0;
              }

              //----------------------
              // Dry deposition (gas)
              //----------------------

              // FIXME: C++ port in progress!
              // mam4::drydep::drydep_xactive(...);

              mam_coupling::vmr2mmr(vmr, adv_mass_kg_per_moles,  q);
              mam_coupling::vmr2mmr(vmrcw, adv_mass_kg_per_moles, qqcw);

              for (int i = offset_aerosol; i < pcnst; ++i) {
                state_q[i] = q[i-offset_aerosol];
                qqcw_long[i] = qqcw[i-offset_aerosol];
              }
              mam4::utils::inject_stateq_to_prognostics(state_q, progs, k);
              mam4::utils::inject_qqcw_to_prognostics(qqcw_long, progs, k);
            });
      });

  // postprocess output
  Kokkos::parallel_for("postprocess", policy, postprocess_);
  Kokkos::fence();
}

void MAMMicrophysics::finalize_impl() {}

}  // namespace scream
