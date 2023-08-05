#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <array>

#include "eamxx_ml_correction_process_interface.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"
#include "share/field/field_utils.hpp"

namespace scream {
// =========================================================================================
MLCorrection::MLCorrection(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  m_ML_model_path = m_params.get<std::string>("ML_model_path");
  m_fields_ml_output_variables = m_params.get<std::vector<std::string>>("ML_output_fields");
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
  constexpr int ps = SCREAM_SMALL_PACK_SIZE;
  m_grid                = grids_manager->get_grid("Physics");
  const auto &grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();  // Number of columns on this rank
  m_num_levs =
      m_grid->get_num_vertical_levels();  // Number of levels per column
  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  FieldLayout scalar3d_layout_mid{{COL, LEV}, {m_num_cols, m_num_levs}};

  /* ----------------------- WARNING --------------------------------*/
  /* The following is a HACK to get things moving, we don't want to
   * add all fields as "updated" long-term.  A separate stream of work
   * is adapting the infrastructure to allow for a generic "add_field" call
   * to be used here which we can then setup using the m_fields_ml_output_variables variable
   */
  add_field<Updated>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  add_field<Updated>("qv",    scalar3d_layout_mid, Q, grid_name, "tracers", ps);
  add_field<Updated>("u",     scalar3d_layout_mid, m/s, grid_name, ps);
  add_field<Updated>("v",     scalar3d_layout_mid, m/s, grid_name, ps);
  /* ----------------------- WARNING --------------------------------*/
  add_group<Updated>("tracers", grid_name, 1, Bundling::Required);
}

void MLCorrection::initialize_impl(const RunType /* run_type */) {
  // Do nothing
}

// =========================================================================================
void MLCorrection::run_impl(const double dt) {
  namespace py      = pybind11;
  // use model time to infer solar zenith angle for the ML prediction
  auto current_ts = timestamp();
  std::string datetime_str = current_ts.get_date_string() + " " + current_ts.get_time_string();
  const auto &qv_field = get_field_out("qv");
  const auto &qv       = qv_field.get_view<Real **, Host>();
  const auto &T_mid_field = get_field_out("T_mid");
  const auto &T_mid       = T_mid_field.get_view<Real **, Host>();
  const auto &u_field = get_field_out("u");
  const auto &u       = u_field.get_view<Real **, Host>();
  const auto &v_field = get_field_out("v");
  const auto &v       = v_field.get_view<Real **, Host>();

  // having m_lat and m_lon in set_grids breaks the standalone unit test
  m_lat  = m_grid->get_geometry_data("lat");
  m_lon  = m_grid->get_geometry_data("lon");
  auto h_lat  = m_lat.get_view<const Real*,Host>();
  auto h_lon  = m_lon.get_view<const Real*,Host>();

  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;
  Int num_tracers = tracers_info->size();
  Real qv_max_before = field_max<Real>(qv_field);
  Real qv_min_before = field_min<Real>(qv_field);
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();  // required for importing numpy
  if ( Py_IsInitialized() == 0 ) {
    py::initialize_interpreter();
  }
  py::module sys = pybind11::module::import("sys");
  sys.attr("path").attr("insert")(1, ML_CORRECTION_CUSTOM_PATH);
  auto py_correction = py::module::import("ml_correction");

  // for qv, we need to stride across number of tracers
  py::object ob1     = py_correction.attr("update_fields")(
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          m_num_cols * m_num_levs, T_mid.data(), py::str{}),
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          m_num_cols * m_num_levs * num_tracers, qv.data(), py::str{}),          
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          m_num_cols * m_num_levs, u.data(), py::str{}),        
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          m_num_cols * m_num_levs, v.data(), py::str{}),       
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          m_num_cols, h_lat.data(), py::str{}),       
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          m_num_cols, h_lon.data(), py::str{}),                                                
      m_num_cols, m_num_levs, num_tracers, dt, m_ML_model_path, datetime_str);
  py::gil_scoped_release no_gil;  
  ekat::enable_fpes(fpe_mask);
  Real qv_max_after = field_max<Real>(qv_field);
  Real qv_min_after = field_min<Real>(qv_field);        
}

// =========================================================================================
void MLCorrection::finalize_impl() {
  // Do nothing
}
// =========================================================================================

}  // namespace scream 
