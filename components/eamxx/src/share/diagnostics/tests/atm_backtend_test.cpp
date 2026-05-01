#include <catch2/catch.hpp>

// Test that the _atm_backtend built-in alias produces numerically identical
// results to the former standalone AtmBackTendDiag.  The alias expands to
// X_minus_X_prev_over_dt, which chains FieldPrevDiag → BinaryOpsDiag → FieldOverDtDiag.

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/diagnostics/field_prev.hpp"
#include "share/diagnostics/field_over_dt.hpp"
#include "share/diagnostics/binary_ops.hpp"
#include "share/field/field_utils.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols,
                                        const int nlevs) {
  const int num_global_cols = ncols * comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names", vos_t{"point_grid"});
  auto &pl = gm_params.sublist("point_grid");
  pl.set<std::string>("type", "point_grid");
  pl.set("aliases", vos_t{"physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm, gm_params);
  gm->build_grids();

  return gm;
}

TEST_CASE("atm_backtend") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);

  util::TimeStamp ts({2024, 1, 1}, {0, 0, 0});

  constexpr int nlevs = 7;
  const int ngcols    = 2 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("physics");

  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());

  Field qc(qc_fid);
  qc.allocate_view();

  int seed = get_random_test_seed(&comm);

  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  // Set time for qc and randomize its values
  qc.get_header().get_tracking().update_time_stamp(ts);
  randomize_uniform(qc, seed++, 0, 200);
  qc.sync_to_dev();

  // Build the composable chain: qc_atm_backtend ≡ qc_minus_qc_prev_over_dt
  // Layer 1: FieldPrevDiag(qc) → output field named "qc_prev"
  ekat::ParameterList prev_params;
  prev_params.set("grid_name", grid->name());
  prev_params.set<std::string>("field_name", "qc");
  auto prev_diag = diag_factory.create("FieldPrevDiag", comm, prev_params);
  prev_diag->set_grids(gm);
  prev_diag->set_required_field(qc);
  prev_diag->initialize(ts, RunType::Initial);

  // The output of FieldPrevDiag is "qc_prev"
  auto qc_prev_field = prev_diag->get_diagnostic();

  // Layer 2: BinaryOpsDiag(qc, minus, qc_prev) → output field "qc_minus_qc_prev"
  ekat::ParameterList minus_params;
  minus_params.set("grid_name", grid->name());
  minus_params.set<std::string>("arg1", "qc");
  minus_params.set<std::string>("arg2", "qc_prev");
  minus_params.set<std::string>("binary_op", "minus");
  auto minus_diag = diag_factory.create("BinaryOpsDiag", comm, minus_params);
  minus_diag->set_grids(gm);
  minus_diag->set_required_field(qc);
  minus_diag->set_required_field(qc_prev_field);
  minus_diag->initialize(ts, RunType::Initial);

  auto qc_minus_qc_prev_field = minus_diag->get_diagnostic();

  // Layer 3: FieldOverDtDiag(qc_minus_qc_prev) → output field "qc_minus_qc_prev_over_dt"
  ekat::ParameterList over_dt_params;
  over_dt_params.set("grid_name", grid->name());
  over_dt_params.set<std::string>("field_name", "qc_minus_qc_prev");
  auto over_dt_diag = diag_factory.create("FieldOverDtDiag", comm, over_dt_params);
  over_dt_diag->set_grids(gm);
  over_dt_diag->set_required_field(qc_minus_qc_prev_field);
  over_dt_diag->initialize(ts, RunType::Initial);

  // First evaluation (before any init_timestep): should throw
  REQUIRE_THROWS(over_dt_diag->compute_diagnostic());

  auto result = over_dt_diag->get_diagnostic();

  const Real a_day = 24.0 * 60.0 * 60.0;  // seconds

  constexpr int ntests = 10;
  constexpr auto tol = std::numeric_limits<Real>::epsilon()*100;

  qc.get_header().get_tracking().update_time_stamp(ts);
  for (int itest = 2; itest < ntests; itest++) {
    // Save current qc before advancing
    auto qc_old = qc.clone();
    qc_old.deep_copy(qc);
    qc_old.sync_to_host();

    // init_timestep: FieldPrevDiag stores current qc; FieldOverDtDiag saves start ts
    prev_diag->init_timestep(ts);
    over_dt_diag->init_timestep(ts);

    // Advance time by one day and update qc
    ts += 86400;
    qc.get_header().get_tracking().update_time_stamp(ts);
    randomize_uniform(qc, seed++, 0, 200);
    qc.sync_to_dev();

    // Evaluate chain in dependency order
    prev_diag->compute_diagnostic();
    minus_diag->compute_diagnostic();
    over_dt_diag->compute_diagnostic();
    result.sync_to_host();

    // Expected: (qc_new - qc_old) / a_day
    auto qc_new_v = qc.get_view<Real**, Host>();
    auto qc_old_v = qc_old.get_view<Real**, Host>();
    auto res_v    = result.get_view<Real**, Host>();

    auto manual = qc.clone();
    manual.update(qc_old,-1,1);
    manual.scale(1/a_day);

    auto diff = result.clone();
    diff.update(manual,-1,1);
    REQUIRE_THAT (inf_norm(diff).as<Real>(), Catch::WithinAbs(0,tol));
  }
}

}  // namespace scream
