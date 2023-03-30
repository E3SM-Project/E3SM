#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <catch2/catch.hpp>
#include <iomanip>

#include "control/atmosphere_driver.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"
#include "physics/ml_correction/atmosphere_ml_correction.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"

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

  EKAT_ASSERT_MSG(dt > 0, "Error! Time step must be positive.\n");

  ekat::Comm atm_comm(MPI_COMM_WORLD);

  auto &proc_factory = AtmosphereProcessFactory::instance();
  auto &gm_factory   = GridsManagerFactory::instance();
  proc_factory.register_product("MLCorrection",
                                &create_atmosphere_process<MLCorrection>);
  gm_factory.register_product("Mesh Free", &create_mesh_free_grids_manager);

  AtmosphereDriver ad;

  ad.initialize(atm_comm, ad_params, t0);

  const auto &grid      = ad.get_grids_manager()->get_grid("Point Grid");
  const auto &field_mgr = *ad.get_field_mgr(grid->name());

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
  Real reference = qv(1, 10);
  reference += 0.1;
  Real reference2 = qv(1, 30);
  reference2 += 0.1;
  pybind11::scoped_interpreter guard{};
  pybind11::module sys = pybind11::module::import("sys");
  sys.attr("path").attr("insert")(1, CUSTOM_SYS_PATH);
  auto py_correction = pybind11::module::import("test_correction");
  py::object ob1     = py_correction.attr("modify_view")(
      py::array_t<Real, py::array::c_style | py::array::forcecast>(
          num_cols * num_levs, qv.data(), py::str{}),
      num_cols, num_levs);
  py::gil_scoped_release no_gil;
  REQUIRE(qv(1, 10) == reference);   // This is the one that is modified
  REQUIRE(qv(1, 30) != reference2);  // This one should be unchanged
  ad.finalize();
}
}  // namespace scream
