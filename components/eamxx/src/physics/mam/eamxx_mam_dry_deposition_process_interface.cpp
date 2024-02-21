#include "physics/mam/eamxx_mam_dry_deposition_process_interface.hpp"

/*
Future work:
Wirte comments
write in/outs for all variables clearly
*/

namespace scream {

// =========================================================================================
MAMDryDep::MAMDryDep(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// =========================================================================================
void MAMDryDep::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto q_unit = kg / kg;
  q_unit.set_string("kg/kg");

  auto n_unit = 1 / kg;  // units of number mixing ratios of tracers
  n_unit.set_string("#/kg");

  auto m3 = m * m * m;  // meter cubed

  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_layout_mid{{COL, LEV}, {ncol_, nlev_}};
  const FieldLayout scalar3d_layout_int{{COL, ILEV}, {ncol_, nlev_ + 1}};

  // Layout for 2D (2d horiz) variable 
  const FieldLayout scalar2d_layout{{COL}, {ncol_}};
  
  // Layout for 4D (2d horiz X 1d vertical x number of modes) variables
  const int num_aero_modes = mam_coupling::num_aero_modes();
  FieldLayout scalar4d_layout_mid{
      {COL, LEV, NMODES}, {ncol_, nlev_, num_aero_modes}};  // mid points

  // -------------------------------------------------------------------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // -------------------------------------------------------------------------------------------------------------------------
  add_field<Required>("T_mid", scalar3d_layout_mid, K,
                      grid_name);  // temperature [K]
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa,
                      grid_name);  // pressure at mid points in [Pa]
  add_field<Required>("p_int", scalar3d_layout_int, Pa,
                      grid_name);  // total pressure
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,
                      grid_name);  // pseudo density in [Pa]
  add_field<Required>("qv", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers");  // specific humidity
  add_field<Required>("qc", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers");  // liquid cloud water [kg/kg] wet
  add_field<Required>("qi", scalar3d_layout_mid, q_unit, grid_name,
                      "tracers");  // ice cloud water [kg/kg] wet
  add_field<Updated>("nc", scalar3d_layout_mid, n_unit, grid_name,
                     "tracers");  // cloud liquid wet number mixing ratio
  add_field<Required>("ni", scalar3d_layout_mid, n_unit, grid_name,
                      "tracers");  // ice number mixing ratio
  add_field<Required>(
      "omega", scalar3d_layout_mid, Pa / s,
      grid_name);  // Vertical pressure velocity [Pa/s] at midpoints

  add_field<Required>("dgncur_awet", scalar4d_layout_mid, m, grid_name);
  add_field<Required>("wetdens", scalar4d_layout_mid, kg / m3, grid_name);
  add_field<Required>("obklen", scalar2d_layout, m, grid_name);
  add_field<Required>("surfric", scalar2d_layout, m / s, grid_name);

  auto nondim = ekat::units::Units::nondimensional();
  add_field<Required>("landfrac", scalar2d_layout, nondim, grid_name);
  add_field<Required>("icefrac", scalar2d_layout, nondim, grid_name);
  add_field<Required>("ocnfrac", scalar2d_layout, nondim, grid_name);
  add_field<Required>("fv", scalar2d_layout, m / s, grid_name);
  add_field<Required>("ram1", scalar2d_layout, s / m, grid_name);

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);

    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name, "tracers");
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, q_unit,
                           grid_name, "tracers");
      }
    }
  }
  // (cloud) aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char *cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);
    // printf("%s \n", int_nmr_field_name);

    add_field<Updated>(cld_nmr_field_name, scalar3d_layout_mid, n_unit,
                       grid_name);
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);

      if(strlen(cld_mmr_field_name) > 0) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_layout_mid, q_unit,
                           grid_name);
      }
    }
  }

  // aerosol-related gases: mass mixing ratios
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_layout_mid, q_unit,
                       grid_name, "tracers");
  }
}

// =========================================================================================
// ON HOST, returns the number of bytes of device memory needed by the above
// Buffer type given the number of columns and vertical levels
size_t MAMDryDep::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

// =========================================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to store
// intermediate (dry) quantities on the given number of columns with the given
// number of vertical levels. Returns the number of bytes allocated.
void MAMDryDep::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMDryDep.");
}

// =========================================================================================
void MAMDryDep::initialize_impl(const RunType run_type) {
  // Gather runtime options
  //(e.g.) runtime_options.lambda_low    = m_params.get<double>("lambda_low");

  wet_atm_.qv    = get_field_in("qv").get_view<const Real **>();
  wet_atm_.qc    = get_field_in("qc").get_view<const Real **>();
  wet_atm_.nc    = get_field_in("nc").get_view<const Real **>();
  wet_atm_.qi    = get_field_in("qi").get_view<const Real **>();
  wet_atm_.ni    = get_field_in("ni").get_view<const Real **>();
  wet_atm_.omega = get_field_in("omega").get_view<const Real **>();

  dry_atm_.T_mid = get_field_in("T_mid").get_view<const Real **>();
  dry_atm_.p_mid = get_field_in("p_mid").get_view<const Real **>();
  dry_atm_.p_del =
      get_field_in("pseudo_density").get_view<const Real **>();
  dry_atm_.qv = buffer_.qv_dry;
  dry_atm_.qc = buffer_.qc_dry;
  dry_atm_.nc = buffer_.nc_dry;
  dry_atm_.qi = buffer_.qi_dry;
  dry_atm_.ni = buffer_.ni_dry;

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
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] =
        get_field_out(gas_mmr_field_name).get_view<Real **>();
    dry_aero_.gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

  // set up our preprocess functor
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_,
                         dry_atm_, dry_aero_);
}

// =========================================================================================
void MAMDryDep::run_impl(const double dt) {

const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();



  /* 
  
  Rough notes:
  
  tair == T_mid
  pmid == p_mid
  pint ==p_int
  pdel == p_del

  state_q: It can be obtained using dry_aero_ and dry_atm_, there is an example of this in the optics code
  qqcw can also be obtained from wet_aero.There is an example in optics code
  Feel free to improvise input variables by adding them to the input.yaml file.



  Function to call in mam4xx/drydep.hpp:
  void aero_model_drydep(
    const ThreadTeam &team,
    const Real fraction_landuse[DryDeposition::n_land_type] ,
    const haero::ConstColumnView tair, const haero::ConstColumnView pmid,
    const haero::ConstColumnView pint, const haero::ConstColumnView pdel,
    const Diagnostics::ColumnTracerView state_q,
    const ColumnView dgncur_awet[AeroConfig::num_modes()],
    const ColumnView wetdens[AeroConfig::num_modes()],
    const Kokkos::View<Real *> qqcw[aero_model::pcnst], const Real obklen,
    const Real ustar, const Real landfrac, const Real icefrac,
    const Real ocnfrac, const Real fricvelin, const Real ram1in,
    const Diagnostics::ColumnTracerView ptend_q,
    bool ptend_lq[aero_model::pcnst], const Real dt,
    const ColumnView aerdepdrycw, const ColumnView aerdepdryis,

    const ColumnView rho,
    const Kokkos::View<Real *> vlc_dry[AeroConfig::num_modes()]
                                      [DryDeposition::aerosol_categories],
    const Kokkos::View<Real *> vlc_trb[AeroConfig::num_modes()]
                                      [DryDeposition::aerosol_categories],
    const Kokkos::View<Real *> vlc_grv[AeroConfig::num_modes()]
                                      [DryDeposition::aerosol_categories],
    const Kokkos::View<Real *> dqdt_tmp[aero_model::pcnst])
  
  */



  std::cout << "End of derydep run" << std::endl;


  
}

// =========================================================================================
}  // namespace scream
