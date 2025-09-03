#include <catch2/catch.hpp>

#include "control/atmosphere_driver.hpp"
#include "physics/register_physics.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

#include <ekat_yaml.hpp>
#include <ekat_fpe.hpp>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <iomanip>

namespace scream {
TEST_CASE("ml_correction-stand-alone", "") {
  using namespace scream;
  using namespace scream::control;
  namespace py      = pybind11;

  std::string fname = "input.yaml";
  ekat::ParameterList ad_params("Atmosphere Driver");
  parse_yaml_file(fname, ad_params);

  const auto& ts     = ad_params.sublist("time_stepping");
  const auto  dt     = ts.get<int>("time_step");
  const auto  nsteps = ts.get<int>("number_of_steps");
  const auto  t0_str = ts.get<std::string>("run_t0");
  const auto  t0     = util::str_to_time_stamp(t0_str);
  const auto  ml     = ad_params.sublist("eamxx").sublist("MLCorrection");
  const auto  ML_model_tq_path = ml.get<std::string>("ml_model_path_tq");
  const auto  ML_model_uv_path = ml.get<std::string>("ml_model_path_uv");

  EKAT_ASSERT_MSG(dt > 0, "Error! Time step must be positive.\n");

  ekat::Comm atm_comm(MPI_COMM_WORLD);

  register_physics();
  register_mesh_free_grids_manager();

  AtmosphereDriver ad;

  ad.initialize(atm_comm, ad_params, t0);

  const auto& grid = ad.get_grids_manager()->get_grid("physics");
  const auto& field_mgr = *ad.get_field_mgr();

  int num_cols = grid->get_num_local_dofs();
  int num_levs = grid->get_num_vertical_levels();

  const auto &qv_field = field_mgr.get_field("qv");
  const auto &qv       = qv_field.get_view<Real **, Host>();

  for(int icol = 0; icol < num_cols; ++icol) {
    for(int jlev = 0; jlev < num_levs; ++jlev) {
      Real phase     = icol * 3.14 / 2.0 / num_cols;
      Real xval      = jlev * 3.14 / 2.0 / num_levs;
      qv(icol, jlev) = (1.0 + std::sin(xval - phase)) / 2.0;
    }
  }
  qv_field.sync_to_dev();
  Real reference = 1e-4;
  int fpe_mask = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();  // required for importing numpy
  if ( Py_IsInitialized() == 0 ) {
    py::initialize_interpreter();
  }
  py::module sys = pybind11::module::import("sys");
  sys.attr("path").attr("insert")(1, CUSTOM_SYS_PATH);
  auto py_correction = py::module::import("test_correction");
  py::object ML_model_tq = py_correction.attr("get_ML_model")(ML_model_tq_path);
  py::object ML_model_uv = py_correction.attr("get_ML_model")(ML_model_uv_path);
  py::object ob1  = py_correction.attr("modify_view")(
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          num_cols * num_levs, qv.data(), py::str{}),
      num_cols, num_levs, ML_model_tq, ML_model_uv);
  py::gil_scoped_release no_gil;
  ekat::enable_fpes(fpe_mask);
  REQUIRE(qv(1, 10) == reference);   // This is the one that is modified
  REQUIRE(qv(0, 10) != reference);  // This one should be unchanged
  ad.finalize();
}

}  // namespace scream
