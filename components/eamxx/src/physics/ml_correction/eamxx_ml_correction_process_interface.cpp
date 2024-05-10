#include "eamxx_ml_correction_process_interface.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"
#include "share/field/field_utils.hpp"

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
  auto Q = kg / kg;
  Q.set_string("kg/kg");
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
  if (not m_ML_correction_unit_test) {
    const auto m2 = m*m;
    const auto s2 = s*s;
    auto Wm2 = W / m / m;
    auto nondim = m/m;
    add_field<Required>("phis", scalar2d, m2/s2, grid_name);
    add_field<Updated>("SW_flux_dn", scalar3d_int, Wm2, grid_name, ps);
    add_field<Required>("sfc_alb_dif_vis", scalar2d, nondim, grid_name);
    add_field<Updated>("sfc_flux_sw_net", scalar2d, Wm2, grid_name);
    add_field<Updated>("sfc_flux_lw_dn", scalar2d, Wm2, grid_name);
    m_lat  = m_grid->get_geometry_data("lat");
    m_lon  = m_grid->get_geometry_data("lon");      
  }

  /* ----------------------- WARNING --------------------------------*/
  /* The following is a HACK to get things moving, we don't want to
   * add all fields as "updated" long-term.  A separate stream of work
   * is adapting the infrastructure to allow for a generic "add_field" call
   * to be used here which we can then setup using the m_fields_ml_output_variables variable
   */
  add_field<Updated>("T_mid",       scalar3d_mid, K,   grid_name, ps);
  add_field<Updated>("qv",          scalar3d_mid, Q,   grid_name, "tracers", ps);
  add_field<Updated>("horiz_winds", vector3d_mid, m/s, grid_name, ps);
  /* ----------------------- WARNING --------------------------------*/
  add_group<Updated>("tracers", grid_name, 1, Bundling::Required);
}

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
}

// =========================================================================================
void MLCorrection::run_impl(const double dt) {
  // use model time to infer solar zenith angle for the ML prediction
  auto current_ts = timestamp();
  std::string datetime_str = current_ts.get_date_string() + " " + current_ts.get_time_string();
  const auto &qv              = get_field_out("qv").get_view<Real **, Host>();
  const auto &T_mid           = get_field_out("T_mid").get_view<Real **, Host>();
  const auto &phis            = get_field_in("phis").get_view<const Real *, Host>();
  const auto &SW_flux_dn      = get_field_out("SW_flux_dn").get_view<Real **, Host>();
  const auto &sfc_alb_dif_vis = get_field_in("sfc_alb_dif_vis").get_view<const Real *, Host>();  
  const auto &sfc_flux_sw_net = get_field_out("sfc_flux_sw_net").get_view<Real *, Host>();
  const auto &sfc_flux_lw_dn  = get_field_out("sfc_flux_lw_dn").get_view<Real *, Host>();
  const auto &u               = get_field_out("horiz_winds").get_component(0).get_view<Real **, Host>();
  const auto &v               = get_field_out("horiz_winds").get_component(1).get_view<Real **, Host>();

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
}

// =========================================================================================
void MLCorrection::finalize_impl() {
  // Do nothing
}
// =========================================================================================

}  // namespace scream 
