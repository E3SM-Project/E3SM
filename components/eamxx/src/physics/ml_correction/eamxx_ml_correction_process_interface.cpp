#include "eamxx_ml_correction_process_interface.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"
#include "share/field/field_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

namespace scream {
// =========================================================================================
MLCorrection::MLCorrection(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  m_ML_model_path_tq = m_params.get<std::string>("ML_model_path_tq");
  m_ML_model_path_uv = m_params.get<std::string>("ML_model_path_uv");
  m_ML_model_path_sfc_fluxes = m_params.get<std::string>("ML_model_path_sfc_fluxes");
  m_fields_ml_output_variables = m_params.get<std::vector<std::string>>("ML_output_fields");
  m_ML_correction_unit_test = m_params.get<bool>("ML_correction_unit_test");
}

// =========================================================================================
void MLCorrection::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  constexpr int ps = Pack::n;
  m_grid                = grids_manager->get_grid("Physics");
  const auto &grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();  // Number of columns on this rank
  m_num_levs =
      m_grid->get_num_vertical_levels();  // Number of levels per column
  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);
  const auto m2 = pow(m,2);
  if (not m_ML_correction_unit_test) {
    const auto s2 = pow(s,2);
    auto nondim = Units::nondimensional();
    add_field<Required>("phis", scalar2d, m2/s2, grid_name);
    add_field<Required>("sfc_alb_dif_vis", scalar2d, nondim, grid_name);
    add_field<Updated>("SW_flux_dn", scalar3d_int, W/m2, grid_name, ps);
    add_field<Updated>("sfc_flux_sw_net", scalar2d, W/m2, grid_name);
    add_field<Updated>("sfc_flux_lw_dn", scalar2d, W/m2, grid_name);
    m_lat  = m_grid->get_geometry_data("lat");
    m_lon  = m_grid->get_geometry_data("lon");
  }

  /* ----------------------- WARNING --------------------------------*/
  /* The following is a HACK to get things moving, we don't want to
   * add all fields as "updated" long-term.  A separate stream of work
   * is adapting the infrastructure to allow for a generic "add_field" call
   * to be used here which we can then setup using the m_fields_ml_output_variables variable
   */
  add_field<Updated>("T_mid",       scalar3d_mid, K,     grid_name, ps);
  add_field<Updated>("horiz_winds", vector3d_mid, m/s,   grid_name, ps);
  /* Note: we also need to update the precipitation after ML commits any column drying */
  add_field<Required>("pseudo_density",      scalar3d_mid, Pa,     grid_name, ps);
  add_field<Updated>("precip_liq_surf_mass", scalar2d,     kg/m2,  grid_name);
  add_field<Updated>("precip_ice_surf_mass", scalar2d,     kg/m2,  grid_name);
  /* ----------------------- WARNING --------------------------------*/
  add_tracer<Updated>("qv", m_grid, kg/kg, ps);
  add_group<Updated>("tracers", grid_name, 1, MonolithicAlloc::Required);
}

// =========================================================================================
void MLCorrection::initialize_impl(const RunType /* run_type */) {
  fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();  // required for importing numpy
  if ( Py_IsInitialized() == 0 ) {
    pybind11::initialize_interpreter();
  }
  pybind11::module sys = pybind11::module::import("sys");
  sys.attr("path").attr("insert")(1, ML_CORRECTION_CUSTOM_PATH);
  py_correction = pybind11::module::import("ml_correction");
  ML_model_tq = py_correction.attr("get_ML_model")(m_ML_model_path_tq);
  ML_model_uv = py_correction.attr("get_ML_model")(m_ML_model_path_uv);
  ML_model_sfc_fluxes = py_correction.attr("get_ML_model")(m_ML_model_path_sfc_fluxes);
  ekat::enable_fpes(fpe_mask);

  // Enforce bounds on quantities adjusted by ML using Field Property Checks
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<LowerBound>(get_field_out("precip_liq_surf_mass"),m_grid,0,true);
  add_postcondition_check<LowerBound>(get_field_out("precip_ice_surf_mass"),m_grid,0,true);
  // Add sanity property check for other variables.
  using Interval = FieldWithinIntervalCheck;
  add_postcondition_check<Interval>(get_field_out("qv"),m_grid,0,0.2,true);
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"),m_grid,-400.0,400.0,false);
}

// =========================================================================================
void MLCorrection::run_impl(const double dt) {
  // use model time to infer solar zenith angle for the ML prediction
  auto current_ts = start_of_step_ts();
  std::string datetime_str = current_ts.get_date_string() + " " + current_ts.get_time_string();

  const auto &phis            = get_field_in("phis").get_view<const Real *, Host>();
  const auto &sfc_alb_dif_vis = get_field_in("sfc_alb_dif_vis").get_view<const Real *, Host>();

  const auto &qv              = get_field_out("qv").get_view<Real **, Host>();
  const auto &T_mid           = get_field_out("T_mid").get_view<Real **, Host>();
  const auto &SW_flux_dn      = get_field_out("SW_flux_dn").get_view<Real **, Host>();
  const auto &sfc_flux_sw_net = get_field_out("sfc_flux_sw_net").get_view<Real *, Host>();
  const auto &sfc_flux_lw_dn  = get_field_out("sfc_flux_lw_dn").get_view<Real *, Host>();
  const auto &u               = get_field_out("horiz_winds").get_component(0).get_view<Real **, Host>();
  const auto &v               = get_field_out("horiz_winds").get_component(1).get_view<Real **, Host>();

  // For precipitation adjustment we need to track the change in column integrated 'qv'
  // So we clone the original qv before ML changes the state so we can back out a qv_tend
  // to use with precip adjustment.
  auto qv_src = get_field_in("qv");
  auto qv_in = qv_src.clone();

  auto h_lat  = m_lat.get_view<const Real*,Host>();
  auto h_lon  = m_lon.get_view<const Real*,Host>();

  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;
  Int num_tracers = tracers_info->size();

  ekat::disable_all_fpes();  // required for importing numpy
  if ( Py_IsInitialized() == 0 ) {
    pybind11::initialize_interpreter();
  }
  // for qv, we need to stride across number of tracers
  pybind11::object ob1     = py_correction.attr("update_fields")(
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols * m_num_levs, T_mid.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols * m_num_levs * num_tracers, qv.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols * m_num_levs, u.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols * m_num_levs, v.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols, h_lat.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols, h_lon.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols, phis.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols * (m_num_levs+1), SW_flux_dn.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols, sfc_alb_dif_vis.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols, sfc_flux_sw_net.data(), pybind11::str{}),
      pybind11::array_t<Real, pybind11::array::c_style | pybind11::array::forcecast>(
          m_num_cols, sfc_flux_lw_dn.data(), pybind11::str{}),
      m_num_cols, m_num_levs, num_tracers, dt,
      ML_model_tq, ML_model_uv, ML_model_sfc_fluxes, datetime_str);
  pybind11::gil_scoped_release no_gil;
  ekat::enable_fpes(fpe_mask);

  // Now back out the qv change abd apply it to precipitation, only if Tq ML is turned on
  if (m_ML_model_path_tq != "None") {
    using PC  = scream::physics::Constants<Real>;
    using KT  = KokkosTypes<DefaultDevice>;
    using MT  = typename KT::MemberType;
    using ESU = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
    const auto &pseudo_density       = get_field_in("pseudo_density").get_view<const Real**>();
    const auto &precip_liq_surf_mass = get_field_out("precip_liq_surf_mass").get_view<Real *>();
    const auto &precip_ice_surf_mass = get_field_out("precip_ice_surf_mass").get_view<Real *>();
    constexpr Real g = PC::gravit;
    const auto num_levs = m_num_levs;
    const auto policy = ESU::get_default_team_policy(m_num_cols, m_num_levs);

    const auto &qv_told = qv_in.get_view<const Real **>();
    const auto &qv_tnew = get_field_in("qv").get_view<const Real **>();
    Kokkos::parallel_for("Compute WVP diff", policy,
                         KOKKOS_LAMBDA(const MT& team) {
      const int icol = team.league_rank();
      auto qold_icol = ekat::subview(qv_told,icol);
      auto qnew_icol = ekat::subview(qv_tnew,icol);
      auto rho_icol  = ekat::subview(pseudo_density,icol);
      Real net_column_moistening = 0;
      // Compute WaterVaporPath Difference
      // The water vapor path (WVP) is calculated as the integral of d_qv over the vertical column
      // which is converted to the units of precipitation which are kg/m*m
      //     WVP = sum( d_qv * pseudo_density / gravity ),
      //       where d_qv = qv_new - qv_old
      //       units sanity check
      //       	d_qv = kg/kg
      //       	pseudo_density = Pa = kg/m/s2
      //       	gravity = m/s2
      //       d_qv * pseduo_density / gravity = kg/kg * kg/m/s2 * s2/m = kg/m2
      // Compute WaterVaporPath Difference
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, num_levs),
                              [&] (const int& ilev, Real& lsum) {
        lsum += (qnew_icol(ilev)-qold_icol(ilev)) * rho_icol(ilev) / g;
      },Kokkos::Sum<Real>(net_column_moistening));
      team.team_barrier();
      // Adjust Precipitation
      //  - Note, we subtract the water vapor path because positive precip represents
      //    a descrease in qv.
      auto tot_precip = precip_liq_surf_mass(icol)+precip_ice_surf_mass(icol);
      if (tot_precip>0) {
        // adjust precip by weighted avg of both phases
        Kokkos::single(Kokkos::PerTeam(team), [&] {
          auto liq_frac = precip_liq_surf_mass(icol)/tot_precip;
          auto ice_frac = precip_ice_surf_mass(icol)/tot_precip;
          precip_liq_surf_mass(icol) -= liq_frac*net_column_moistening;
          precip_ice_surf_mass(icol) -= ice_frac*net_column_moistening;
	});
      } else {
        // Apply all the adjustment to a single phase based on surface temperature
        Kokkos::single(Kokkos::PerTeam(team), [&] {
          auto T_icol = ekat::subview(T_mid,icol);
          if (T_icol(num_levs-1)>273.15) {
            precip_liq_surf_mass(icol) -= net_column_moistening;
          } else {
            precip_ice_surf_mass(icol) -= net_column_moistening;
          }
	});
      }
      // Note, with the above formulation it is possible for the precipitation to go negative.
      // We rely on the field property checker defined in the intialization function to repair
      // any instances of negative precipitation.
    });
  }
}

// =========================================================================================
void MLCorrection::finalize_impl() {
  // Do nothing
}
// =========================================================================================

}  // namespace scream
