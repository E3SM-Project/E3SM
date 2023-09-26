#include <physics/mam/eamxx_mam_microphysics_process_interface.hpp>
#include <share/property_checks/field_lower_bound_check.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

#include "scream_config.h" // for SCREAM_CIME_BUILD

#include <ekat/ekat_assert.hpp>

// NOTE: see the impl/ directory for the contents of the impl namespace

namespace scream
{

MAMMicrophysics::MAMMicrophysics(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params),
    aero_config_() {
  configure(params);
}

AtmosphereProcessType MAMMicrophysics::type() const {
  return AtmosphereProcessType::Physics;
}

std::string MAMMicrophysics::name() const {
  return "mam4_micro";
}

void MAMMicrophysics::configure(const ekat::ParameterList& params) {
  // enable/disable specific parameterizations
  config_.do_cond = true;
  config_.do_rename = true;
  config_.do_newnuc = true;
  config_.do_coag = true;

  // configure the specific aerosol microphysics parameterizations

  config_.nucleation = {};
  config_.nucleation.dens_so4a_host = 1770.0;
  config_.nucleation.mw_so4a_host = 115.0;
  config_.nucleation.newnuc_method_user_choice = 2;
  config_.nucleation.pbl_nuc_wang2008_user_choice = 1;
  config_.nucleation.adjust_factor_pbl_ratenucl = 1.0;
  config_.nucleation.accom_coef_h2so4 = 1.0;
  config_.nucleation.newnuc_adjust_factor_dnaitdt = 1.0;

  // these parameters guide the coupling between parameterizations
  // NOTE: mam4xx was ported with these parameters set thus. It's probably not
  // NOTE: safe to change these without code modifications.
  config_.gaexch_h2so4_uptake_optaa = 2;
  config_.newnuc_h2so4_conc_optaa = 2;
}

void MAMMicrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  auto q_unit = kg/kg; // mass mixing ratios [kg stuff / kg air]
  q_unit.set_string("kg/kg");
  auto n_unit = 1/kg;  // number mixing ratios [# / kg air]
  n_unit.set_string("#/kg");
  Units nondim(0,0,0,0,0,0,0);

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();      // number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // number of levels per column

  // define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{ {COL}, {ncol_} };

  // layout for 3D (2d horiz X 1d vertical) variables
  FieldLayout scalar3d_layout_mid{ {COL, LEV}, {ncol_, nlev_} };

  // define fields needed in mam4xx

  // atmospheric quantities
  add_field<Required>("omega", scalar3d_layout_mid, Pa/s, grid_name); // vertical pressure velocity
  add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name); // Temperature
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name); // total pressure
  add_field<Required>("qv", scalar3d_layout_mid, q_unit, grid_name, "tracers"); // specific humidity
  add_field<Required>("qi", scalar3d_layout_mid, q_unit, grid_name, "tracers"); // ice wet mixing ratio
  add_field<Required>("ni", scalar3d_layout_mid, n_unit, grid_name, "tracers"); // ice number mixing ratio
  add_field<Required>("pbl_height", scalar2d_layout_col, m, grid_name); // planetary boundary layer height
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name); // p_del, hydrostatic pressure
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name); // cloud fraction

  // droplet activation can alter cloud liquid and number mixing ratios
  add_field<Updated>("qc", scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud liquid wet mixing ratio
  add_field<Updated>("nc", scalar3d_layout_mid, n_unit, grid_name, "tracers"); // cloud liquid wet number mixing ratio

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit, grid_name, "tracers");
    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char* int_mmr_field_name = mam_coupling::int_aero_mmr_field_name(m, a);
      if (strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name, "tracers");
      }
    }
  }

  // aerosol-related gases: mass mixing ratios
  for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char* gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name, "tracers");
  }

  // Tracers group -- do we need this in addition to the tracers above? In any
  // case, this call should be idempotent, so it can't hurt.
  add_group<Updated>("tracers", grid_name, 1, Bundling::Required);
}

// this checks whether we have the tracers we expect
void MAMMicrophysics::
set_computed_group_impl(const FieldGroup& group) {
  const auto& name = group.m_info->m_group_name;
  EKAT_REQUIRE_MSG(name=="tracers",
    "Error! MAM4 expects a 'tracers' field group (got '" << name << "')\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
    "Error! MAM4 expects bundled fields for tracers.\n");

  // how many aerosol/gas tracers do we expect?
  int num_tracers = 2 * (mam_coupling::num_aero_modes() +
                         mam_coupling::num_aero_tracers()) +
                    mam_coupling::num_aero_gases();
  EKAT_REQUIRE_MSG(group.m_info->size() >= num_tracers,
    "Error! MAM4 requires at least " << num_tracers << " aerosol tracers.");
}

size_t MAMMicrophysics::requested_buffer_size_in_bytes() const
{
  return mam_coupling::buffer_size(ncol_, nlev_);
}

void MAMMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Insufficient buffer size.\n");

  size_t used_mem = mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMMicrophysics.");
}

void MAMMicrophysics::initialize_impl(const RunType run_type) {

  step_ = 0;

  // populate the wet and dry atmosphere states with views from fields and
  // the buffer
  wet_atm_.qv = get_field_in("qv").get_view<const Real**>();
  wet_atm_.qc = get_field_out("qc").get_view<Real**>();
  wet_atm_.nc = get_field_out("nc").get_view<Real**>();
  wet_atm_.qi = get_field_in("qi").get_view<const Real**>();
  wet_atm_.ni = get_field_in("ni").get_view<const Real**>();
  wet_atm_.omega = get_field_in("omega").get_view<const Real**>();

  dry_atm_.T_mid     = get_field_in("T_mid").get_view<const Real**>();
  dry_atm_.p_mid     = get_field_in("p_mid").get_view<const Real**>();
  dry_atm_.p_del     = get_field_in("pseudo_density").get_view<const Real**>();
  dry_atm_.cldfrac   = get_field_in("cldfrac_tot").get_view<const Real**>(); // FIXME: tot or liq?
  dry_atm_.pblh      = get_field_in("pbl_height").get_view<const Real*>();
  dry_atm_.z_mid     = buffer_.z_mid;
  dry_atm_.dz        = buffer_.dz;
  dry_atm_.z_iface   = buffer_.z_iface;
  dry_atm_.qv        = buffer_.qv_dry;
  dry_atm_.qc        = buffer_.qc_dry;
  dry_atm_.nc        = buffer_.nc_dry;
  dry_atm_.qi        = buffer_.qi_dry;
  dry_atm_.ni        = buffer_.ni_dry;
  dry_atm_.w_updraft = buffer_.w_updraft;
  dry_atm_.z_surf = 0.0; // FIXME: for now

  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;

  // perform any initialization work
  if (run_type==RunType::Initial) {
  }

  // set wet/dry aerosol state data (interstitial aerosols only)
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    wet_aero_.int_aero_nmr[m] = get_field_out(int_nmr_field_name).get_view<Real**>();
    dry_aero_.int_aero_nmr[m] = buffer_.dry_int_aero_nmr[m];
    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char* int_mmr_field_name = mam_coupling::int_aero_mmr_field_name(m, a);
      if (strlen(int_mmr_field_name) > 0) {
        wet_aero_.int_aero_mmr[m][a] = get_field_out(int_mmr_field_name).get_view<Real**>();
        dry_aero_.int_aero_mmr[m][a] = buffer_.dry_int_aero_mmr[m][a];
      }
    }
  }

  // set wet/dry aerosol-related gas state data
  for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char* mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] = get_field_out(mmr_field_name).get_view<Real**>();
    dry_aero_.gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

  // set up our preprocess/postprocess functors
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_, dry_aero_);
  postprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_, dry_aero_);

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
  // FIXME: we'll probably need this later, but we'll just use ATMBufferManager for now
  //const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  //workspace_mgr_.setup(buffer_.wsm_data, nlev_+1, 13+(n_wind_slots+n_trac_slots), default_policy);
}

void MAMMicrophysics::run_impl(const double dt) {

  const auto scan_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);
  const auto policy      = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  // reset internal WSM variables
  //workspace_mgr_.reset_internals();

  // FIXME: nothing depends on simulation time (yet), so we can just use zero for now
  double t = 0.0;

  // loop over atmosphere columns and compute aerosol microphyscs
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    const Int icol = team.league_rank(); // column index

    // fetch column-specific atmosphere state data
    auto atm = mam_coupling::atmosphere_for_column(dry_atm_, icol);

    // set surface state data
    haero::Surface sfc{};

    // fetch column-specific subviews into aerosol prognostics
    mam4::Prognostics progs = mam_coupling::interstitial_aerosols_for_column(dry_aero_, icol);

    // set up diagnostics
    mam4::Diagnostics diags(nlev_);

    // compute aerosol microphysics on each vertical level within this column
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_), [&](const int k) {

      constexpr int num_modes = mam4::AeroConfig::num_modes();
      constexpr int num_aero_ids = mam4::AeroConfig::num_aerosol_ids();
      constexpr int gas_pcnst = mam_coupling::gas_pcnst();
      constexpr int nqtendbb = mam_coupling::nqtendbb();

      // extract atm state variables (input)
      Real temp = atm.temperature(k);
      Real pmid = atm.pressure(k);
      Real pdel = atm.hydrostatic_dp(k);
      Real zm   = atm.height(k);
      Real pblh = atm.planetary_boundary_layer_height;
      Real qv   = atm.vapor_mixing_ratio(k);
      Real cld  = atm.cloud_fraction(k);

      // aerosol/gas species work arrays (volume mixing ratios (VMR), equivalent
      // to tracer mixing ratios (TMR))
      Real q[gas_pcnst] = {};
      Real qqcw[gas_pcnst] = {};

      // extract aerosol state variables into "working arrays"
      // (in EAM, this is done in the gas_phase_chemdr subroutine defined within
      //  mozart/mo_gas_phase_chemdr.F90)
      mam_coupling::transfer_prognostics_to_work_arrays(progs, k, q, qqcw);

      // aerosol/gas species tendencies (output)
      Real q_tendbb[gas_pcnst][nqtendbb] = {};
      Real qqcw_tendbb[gas_pcnst][nqtendbb] = {};

      //-------------------------------------------------------------------
      // Chemistry (doesn't belong here, but it's here for now to minimize
      // necessary additional tracers like the work array copies below)
      //-------------------------------------------------------------------

      // create work array copies to retain "pre-chemistry" values
      Real q_pregaschem[gas_pcnst] = {};
      Real q_precldchem[gas_pcnst] = {};
      Real qqcw_precldchem[gas_pcnst] = {};
      for (int i = 0; i < gas_pcnst; ++i) {
        q_pregaschem[i] = q[i];
        q_precldchem[i] = q[i];
        qqcw_precldchem[i] = qqcw[i];
      }

      // FIXME: Chemistry goes here

      //----------------------
      // Aerosol microphysics
      //----------------------

      // calculate aerosol water content using water uptake treatment
      // * dry and wet diameters [m]
      // * wet densities [kg/m3]
      // * aerosol water mass mixing ratio [kg/kg]
      Real dgncur_a[num_modes]    = {};
      Real dgncur_awet[num_modes] = {};
      Real wetdens[num_modes]     = {};
      Real qaerwat[num_modes]     = {};
      {
        // get some information about aerosol species
        // FIXME: this isn't great!
        constexpr int maxd_aspectype = mam4::water_uptake::maxd_aspectype;
        int nspec_amode[num_modes], lspectype_amode[maxd_aspectype][num_modes];
        Real specdens_amode[maxd_aspectype], spechygro[maxd_aspectype];
        mam4::water_uptake::get_e3sm_parameters(nspec_amode, lspectype_amode,
          specdens_amode, spechygro);

        // extract aerosol tracers for this level into state_q, which is needed
        // for computing dry aerosol properties below
        // FIXME: we should eliminate this index translation stuff
        constexpr int nvars = mam4::water_uptake::nvars;
        Real state_q[nvars]; // aerosol tracers for level k
        for (int imode = 0; imode < num_modes; ++imode) {
          int la, lc; // interstitial and cloudborne indices within state_q

          // number mixing ratios
          mam4::convproc::assign_la_lc(imode, -1, la, lc);
          state_q[la] = progs.n_mode_i[imode](k);
          state_q[lc] = progs.n_mode_c[imode](k);
          // aerosol mass mixing ratios
          for (int iaero = 0; iaero < num_aero_ids; ++iaero) {
            mam4::convproc::assign_la_lc(imode, iaero, la, lc);
            auto mode = static_cast<mam4::ModeIndex>(imode);
            auto aero = static_cast<mam4::AeroId>(iaero);
            int ispec = mam4::aerosol_index_for_mode(mode, aero);
            if (ispec != -1) {
              state_q[la] = progs.q_aero_i[imode][ispec](k);
              state_q[lc] = progs.q_aero_c[imode][ispec](k);
            }
          }
        }

        // compute the dry volume for each mode, and from it the current dry
        // geometric nominal particle diameter.
        // FIXME: We have to do some gymnastics here to set up the calls to
        // FIXME: calcsize. This could be improved.
        Real inv_densities[num_modes][num_aero_ids] = {};
        for (int imode = 0; imode < num_modes; ++imode) {
          const int n_spec = num_species_mode(imode);
          for (int ispec = 0; ispec < n_spec; ++ispec) {
            const int iaer = static_cast<int>(mode_aero_species(imode, ispec));
            const Real density = aero_species(iaer).density;
            inv_densities[imode][ispec] = 1.0 / density;
          }
        }
        for (int imode = 0; imode < num_modes; ++imode) {
          Real dryvol_i, dryvol_c; // interstitial and cloudborne dry volumes
          mam4::calcsize::compute_dry_volume_k(k, imode, inv_densities, progs,
            dryvol_i, dryvol_c);

          // NOTE: there's some disagreement over whether vol2num should be called
          // NOTE: num2vol here, so I'm just adopting the nomenclature used by
          // NOTE: the following call to calcsize)
          const mam4::Mode& mode = mam4::modes(imode);
          Real vol2num_min = 1.0/mam4::conversions::mean_particle_volume_from_diameter(
              mode.max_diameter, mode.mean_std_dev);
          Real vol2num_max = 1.0/mam4::conversions::mean_particle_volume_from_diameter(
              mode.min_diameter, mode.mean_std_dev);
          Real vol2num;
          mam4::calcsize::update_diameter_and_vol2num(dryvol_i,
            progs.n_aero_i[imode], vol2num_min, vol2num_max,
            mode.min_diameter, mode.max_diameter, mode.mean_std_dev,
            dgncur_a[imode], vol2num);
        }

        // calculate dry aerosol properties
        Real hygro[num_modes], naer[num_modes], dryrad[num_modes],
             dryvol[num_modes], drymass[num_modes],
             rhcrystal[num_modes], rhdeliques[num_modes], specdens_1[num_modes];
        mam4::water_uptake::modal_aero_wateruptake_dryaer(nspec_amode, specdens_amode,
          spechygro, lspectype_amode, state_q, dgncur_a, hygro,
          naer, dryrad, dryvol, drymass, rhcrystal, rhdeliques, specdens_1);

        // calculate wet aerosol properties
        Real rh = mam4::conversions::relative_humidity_from_vapor_mixing_ratio(qv, temp, pmid);
        Real wetrad[num_modes], wetvol[num_modes], wtrvol[num_modes];
        mam4::water_uptake::modal_aero_wateruptake_wetaer(rhcrystal, rhdeliques, dgncur_a,
          dryrad, hygro, rh, naer, dryvol, wetrad, wetvol, wtrvol, dgncur_awet,
          qaerwat);
        mam4::water_uptake::modal_aero_wateruptake_wetdens(wetvol, wtrvol,
          drymass, specdens_1, wetdens);
      }

      // compute aerosol microphysics
      impl::modal_aero_amicphys_intr(config_, ncol_, step_, dt, t, pmid, pdel,
                                     zm, pblh, qv, cld, q, qqcw, q_pregaschem,
                                     q_precldchem, qqcw_precldchem, q_tendbb,
                                     qqcw_tendbb, dgncur_a, dgncur_awet,
                                     wetdens, qaerwat);

      // unpack updated prognostics from work arrays
      mam_coupling::transfer_work_arrays_to_prognostics(q, qqcw, progs, k);
    });
  });

  // postprocess output
  Kokkos::parallel_for("postprocess", policy, postprocess_);
  Kokkos::fence();
}

void MAMMicrophysics::finalize_impl()
{
}

} // namespace scream
