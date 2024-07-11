#include "physics/mam/eamxx_mam_dry_deposition_process_interface.hpp"

#include "mam4xx/drydep.hpp"

/*
Inputs:

Atmosphere:
  temperature        "T_mid"
  pressure           "p_mid"
  interface_pressure "p_int"
  hydrostatic_dp     "pseudo_density"

Diagnostics:
  tracer_mixing_ratio              "qtracers"
  wet_geometric_mean_diameter_i    "dgncur_awet"
  wet_density                      "wetdens"

Diagnostic Scalar Parameters, one per column:
  "Obukhov_length"
  "surface_friction_velocty"
  "land_fraction"
  "ice_fraction"
  "ocean_fraction"
  "friction_velocity"
  "aerodynamical_resistance"

Prognostics:
  n_mode_c    "mam_coupling::cld_aero_nmr_field_name(m)"
  q_aero_c    "mam_coupling::int_aero_nmr_field_name(m)"



Outputs:

Diagnostics:
  d_tracer_mixing_ratio_dt                 "d_qtracers_dt"
  deposition_flux_of_cloud_borne_aerosols
"deposition_flux_of_cloud_borne_aerosols"
  deposition_flux_of_interstitial_aerosols
"deposition_flux_of_interstitial_aerosols"

Computed internally, could be exposed if needed by another process
  vlc_grv(nlev) : dep velocity of gravitational settling [m/s]
  vlc_trb(nlev) : dep velocity of turbulent dry deposition [m/s]
  vlc_dry(nlev) : dep velocity, sum of vlc_grv and vlc_trb [m/s]

Tendencies:
  n_mode_c   "Tendencies"
  q_aero_c   "Tendencies"

*/

namespace scream {

MAMDryDep::MAMDryDep(const ekat::Comm &comm, const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMDryDep::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  // set grid for all the inputs and outputs
  // use physics grid
  grid_ = grids_manager->get_grid("Physics");

  // Name of the grid
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (2d horiz) variable
  // FIXME: Fix layouts based on new format
  const FieldLayout scalar2d{{COL}, {ncol_}};

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_mid{{COL, LEV}, {ncol_, nlev_}};
  const FieldLayout scalar3d_int{{COL, ILEV}, {ncol_, nlev_ + 1}};

  // Layout for 4D (2d horiz X 1d vertical x number of modes) variables
  // at mid points
  const int num_aero_modes = mam_coupling::num_aero_modes();
  FieldLayout scalar4d_mid{{NMODES, COL, LEV}, {num_aero_modes, ncol_, nlev_}};

  // Layout for tracers.
  const int pcnst = mam4::aero_model::pcnst;
  FieldLayout scalar4d_q{{COL, LEV, CMP}, {ncol_, nlev_, pcnst}};
  FieldLayout scalar4d_qqcw_tends{{COL, CMP, LEV}, {ncol_, pcnst, nlev_}};

  using namespace ekat::units;

  auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  auto n_unit = 1 / kg;   // units of number mixing ratios of tracers

  auto nondim = ekat::units::Units::nondimensional();

  auto m3 = m * m * m;  // meter cubed

  // --------------------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // --------------------------------------------------------------------------

  // ----------- Atmospheric quantities -------------
  // Specific humidity [kg/kg]
  add_field<Required>("qv", scalar3d_mid, q_unit, grid_name, "tracers");

  // Cloud liquid mass mixing ratio [kg/kg]
  add_field<Required>("qc", scalar3d_mid, q_unit, grid_name, "tracers");

  // Cloud ice mass mixing ratio [kg/kg]
  add_field<Required>("qi", scalar3d_mid, q_unit, grid_name, "tracers");

  // Cloud liquid number mixing ratio [1/kg]
  add_field<Required>("nc", scalar3d_mid, n_unit, grid_name, "tracers");

  // Cloud ice number mixing ratio [1/kg]
  add_field<Required>("ni", scalar3d_mid, n_unit, grid_name, "tracers");

  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);

  // Vertical pressure velocity [Pa/s] at midpoints
  add_field<Required>("omega", scalar3d_mid, Pa / s, grid_name);

  // Total pressure [Pa] at midpoints
  add_field<Required>("p_mid", scalar3d_mid, Pa, grid_name);

  // Total pressure [Pa] at interfaces
  add_field<Required>("p_int", scalar3d_int, Pa, grid_name);

  // Layer thickness(pdel) [Pa] at midpoints
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);

  // Planetary boundary layer height [m]
  add_field<Required>("pbl_height", scalar2d, m, grid_name);

  static constexpr auto m2 = m * m;
  static constexpr auto s2 = s * s;

  // Surface geopotential [m2/s2]
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  //----------- Variables from microphysics scheme -------------

  // Total cloud fraction [fraction]
  // FIXME: Is is cldfrac_liq instead? find out
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);

  //----------- Variables from coupler (land component)---------
  // Obukhov length [m]
  add_field<Required>("Obukhov_length", scalar2d, m, grid_name);

  // Surface friction velocty [m]
  add_field<Required>("surface_friction_velocty", scalar2d, m / s, grid_name);

  // Land fraction [fraction]
  add_field<Required>("land_fraction", scalar2d, nondim, grid_name);

  // Friction velocity from land model [m/s]
  add_field<Required>("friction_velocity", scalar2d, m / s, grid_name);

  // Aerodynamical resistance from land model [s/m]
  add_field<Required>("aerodynamical_resistance", scalar2d, s / m, grid_name);

  //----------- Variables from coupler (ice component)---------

  // Ice fraction [unitless]
  add_field<Required>("ice_fraction", scalar2d, nondim, grid_name);

  //----------- Variables from coupler (ocean component)---------
  // Ocean fraction [unitless]
  add_field<Required>("ocean_fraction", scalar2d, nondim, grid_name);

  // ---------------------------------------------------------------------
  // These variables are "updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------

  add_field<Updated>("wetdens", scalar4d_mid, kg / m3, grid_name);
  add_field<Updated>("dgncur_awet", scalar4d_mid, m, grid_name);

  // Tempary memory to format the wet_aero.int_aero_nmr and
  // wet_aero.int_aero_mmr views into the order expected by mam4xx.
  add_field<Computed>("qtracers", scalar4d_q, kg / kg, grid_name);
  add_field<Computed>("d_qtracers_dt", scalar4d_q, kg / kg / s, grid_name);

  // surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]
  add_field<Computed>("deposition_flux_of_cloud_borne_aerosols", scalar3d_mid,
                      1 / (m * m) / s, grid_name);
  // surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
  add_field<Computed>("deposition_flux_of_interstitial_aerosols", scalar3d_mid,
                      1 / (m * m) / s, grid_name);
  add_field<Computed>("Tendencies", scalar4d_qqcw_tends, kg / kg / s,
                      grid_name);

  // TODO: The following are not used by drydep but to create a dry atmosphere
  // object.

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  for(int m = 0; m < num_aero_modes; ++m) {
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);

    add_field<Updated>(int_nmr_field_name, scalar3d_mid, n_unit, grid_name,
                       "tracers");
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_mid, q_unit, grid_name,
                           "tracers");
      }
    }
  }
  // (cloud) aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for(int m = 0; m < num_aero_modes; ++m) {
    const char *cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);

    add_field<Updated>(cld_nmr_field_name, scalar3d_mid, n_unit, grid_name);
    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);

      if(strlen(cld_mmr_field_name) > 0) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_mid, q_unit, grid_name);
      }
    }
  }

  // aerosol-related gases: mass mixing ratios
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_mid, q_unit, grid_name,
                       "tracers");
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

  wet_atm_.qv = get_field_in("qv").get_view<const Real **>();
  wet_atm_.qc = get_field_in("qc").get_view<const Real **>();
  wet_atm_.nc = get_field_in("nc").get_view<const Real **>();
  wet_atm_.qi = get_field_in("qi").get_view<const Real **>();
  wet_atm_.ni = get_field_in("ni").get_view<const Real **>();

  dry_atm_.T_mid   = get_field_in("T_mid").get_view<const Real **>();
  dry_atm_.p_mid   = get_field_in("p_mid").get_view<const Real **>();
  dry_atm_.p_del   = get_field_in("pseudo_density").get_view<const Real **>();
  dry_atm_.cldfrac = get_field_in("cldfrac_tot")
                         .get_view<const Real **>();  // FIXME: tot or liq?
  dry_atm_.pblh = get_field_in("pbl_height").get_view<const Real *>();

  obukhov_length_ = get_field_in("Obukhov_length").get_view<const Real *>();
  surface_friction_velocty_ =
      get_field_in("surface_friction_velocty").get_view<const Real *>();
  land_fraction_  = get_field_in("land_fraction").get_view<const Real *>();
  ice_fraction_   = get_field_in("ice_fraction").get_view<const Real *>();
  ocean_fraction_ = get_field_in("ocean_fraction").get_view<const Real *>();
  friction_velocity_ =
      get_field_in("friction_velocity").get_view<const Real *>();
  aerodynamical_resistance_ =
      get_field_in("aerodynamical_resistance").get_view<const Real *>();

  dry_atm_.p_int     = get_field_in("p_int").get_view<const Real **>();
  dry_atm_.z_mid     = buffer_.z_mid;
  dry_atm_.z_iface   = buffer_.z_iface;
  dry_atm_.dz        = buffer_.dz;
  dry_atm_.qv        = buffer_.qv_dry;
  dry_atm_.qc        = buffer_.qc_dry;
  dry_atm_.nc        = buffer_.nc_dry;
  dry_atm_.qi        = buffer_.qi_dry;
  dry_atm_.ni        = buffer_.ni_dry;
  dry_atm_.w_updraft = buffer_.w_updraft;
  dry_atm_.z_surf    = 0.0;  // FIXME: for now

  dgncur_awet_ = get_field_out("dgncur_awet").get_view<Real ***>();
  wet_dens_    = get_field_out("wetdens").get_view<Real ***>();

  qtracers_      = get_field_out("qtracers").get_view<Real ***>();
  d_qtracers_dt_ = get_field_out("d_qtracers_dt").get_view<Real ***>();

  aerdepdrycw_ = get_field_out("deposition_flux_of_cloud_borne_aerosols")
                     .get_view<Real **>();
  aerdepdryis_ = get_field_out("deposition_flux_of_interstitial_aerosols")
                     .get_view<Real **>();
  tendencies_ = get_field_out("Tendencies").get_view<Real ***>();
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
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                         dry_aero_);
}
namespace {
void compute_tendencies(
    const int ncol, const int nlev, const mam4::DryDeposition dry_deposition,
    const double dt, const MAMDryDep::const_view_1d obklen,
    const MAMDryDep::const_view_1d surfric,
    const MAMDryDep::const_view_1d landfrac,
    const MAMDryDep::const_view_1d icefrac,
    const MAMDryDep::const_view_1d ocnfrac,
    const MAMDryDep::const_view_1d friction_velocity,
    const MAMDryDep::const_view_1d aerodynamical_resistance,
    MAMDryDep::view_3d qtracers, MAMDryDep::view_3d d_qtracers_dt,
    const MAMDryDep::view_3d dgncur_awet_, const MAMDryDep::view_3d wet_dens_,
    const mam_coupling::DryAtmosphere dry_atm,
    const mam_coupling::AerosolState dry_aero,
    const mam_coupling::AerosolState wet_aero, MAMDryDep::view_2d aerdepdrycw,
    MAMDryDep::view_2d aerdepdryis, MAMDryDep::view_3d tendencies) {
  const auto policy =
      ekat::ExeSpaceUtils<MAMDryDep::KT::ExeSpace>::get_default_team_policy(
          ncol, nlev);
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const MAMDryDep::KT::MemberType &team) {
        const int num_aero_modes   = mam_coupling::num_aero_modes();
        const int num_aero_species = mam_coupling::num_aero_species();
        const int icol             = team.league_rank();
        const Real t               = 0;

        compute_wet_mixing_ratios(team, dry_atm, dry_aero, wet_aero, icol);
        team.team_barrier();
        Kokkos::parallel_for(
            Kokkos::TeamVectorRange(team, nlev), [&](const int lev) {
              for(int mode = 0; mode < num_aero_modes; ++mode) {
                int icnst = mam4::ConvProc::numptrcw_amode(mode);
                qtracers(icol, lev, icnst) =
                    wet_aero.int_aero_nmr[mode](icol, lev);
                for(int species = 0; species < num_aero_species; ++species) {
                  icnst = mam4::ConvProc::lmassptrcw_amode(species, mode);
                  if(-1 < icnst) {
                    qtracers(icol, lev, icnst) =
                        wet_aero.int_aero_mmr[mode][species](icol, lev);
                  }
                }
              }
            });
        team.team_barrier();

        mam4::Atmosphere atm    = atmosphere_for_column(dry_atm, icol);
        mam4::Prognostics progs = aerosols_for_column(dry_aero, icol);
        mam4::Surface surf;
        mam4::Diagnostics diags;
        mam4::Tendencies tends;

        for(int i = 0; i < num_aero_modes; ++i) {
          diags.wet_geometric_mean_diameter_i[i] =
              ekat::subview(dgncur_awet_, i, icol);
          diags.wet_density[i] = ekat::subview(wet_dens_, i, icol);
        }
        diags.tracer_mixing_ratio      = ekat::subview(qtracers, icol);
        diags.d_tracer_mixing_ratio_dt = ekat::subview(d_qtracers_dt, icol);
        diags.deposition_flux_of_cloud_borne_aerosols =
            ekat::subview(aerdepdrycw, icol);
        diags.deposition_flux_of_interstitial_aerosols =
            ekat::subview(aerdepdryis, icol);

        diags.Obukhov_length           = obklen[icol];
        diags.surface_friction_velocty = surfric[icol];
        diags.land_fraction            = landfrac[icol];
        diags.ice_fraction             = icefrac[icol];
        diags.ocean_fraction           = ocnfrac[icol];
        diags.friction_velocity        = friction_velocity[icol];
        diags.aerodynamical_resistance = aerodynamical_resistance[icol];

        // Fill Tendency views
        for(int m = 0; m < num_aero_modes; ++m) {
          int iconv         = mam4::ConvProc::numptrcw_amode(m);
          tends.n_mode_c[m] = ekat::subview(tendencies, icol, iconv);
          for(int a = 0; a < num_aero_species; ++a) {
            iconv = mam4::ConvProc::lmassptrcw_amode(a, m);
            if(-1 < iconv)
              tends.q_aero_c[m][a] = ekat::subview(tendencies, icol, iconv);
          }
        }

        const mam4::AeroConfig aero_config;
        dry_deposition.compute_tendencies(aero_config, team, t, dt, atm, surf,
                                          progs, diags, tends);
      });
}
}  // namespace

// =========================================================================================
void MAMDryDep::run_impl(const double dt) {
  using DryDep           = mam4::DryDeposition;
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  const DryDep::Config process_config;
  // Future enhancement:
  // This is where we can set the fraction of land use parameters:
  // if (do_soilw) {
  //   const Real *landuse = get_landuse_from_file(do_soilw);
  //   for (int i=0; i<DryDep::n_land_type; ++i)
  //     process_config.fraction_landuse[i] = landuse[i];
  // }

  DryDep dry_deposition;
  const mam4::AeroConfig aero_config;
  dry_deposition.init(aero_config, process_config);

  // validate is not implimented as yet.
  // const auto policy =
  // ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  // Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
  //   dry_deposition.validate( aero_config, team, atm, progs);
  // }
  // Kokkos::fence();

  compute_tendencies(ncol_, nlev_, dry_deposition, dt, obukhov_length_,
                     surface_friction_velocty_, land_fraction_, ice_fraction_,
                     ocean_fraction_, friction_velocity_,
                     aerodynamical_resistance_, qtracers_,
                     d_qtracers_dt_,  // d_qtracers_dt_ is an output
                     dgncur_awet_, wet_dens_, dry_atm_, dry_aero_, wet_aero_,
                     // Outputs:
                     aerdepdrycw_, aerdepdryis_, tendencies_);
  Kokkos::fence();
}
// =========================================================================================
}  // namespace scream
