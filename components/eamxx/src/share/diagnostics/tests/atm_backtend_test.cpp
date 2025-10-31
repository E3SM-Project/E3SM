#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
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

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 7;
  const int ngcols    = 2 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("physics");

  // Input (randomized) qc
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());

  Field qc(qc_fid);
  qc.allocate_view();

  // Random number generator seed
  int seed = get_random_test_seed(&comm);

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("AtmBackTendDiag", comm,
                                     params));  // No 'tendency_name'

  // Set time for qc and randomize its values
  qc.get_header().get_tracking().update_time_stamp(t0);
  randomize_uniform(qc, seed++, 0, 200);

  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("tendency_name", "qc");
  auto diag = diag_factory.create("AtmBackTendDiag", comm, params);
  diag->set_grids(gm);
  diag->set_required_field(qc);
  diag->initialize(t0, RunType::Initial);

  // Run diag
  diag->compute_diagnostic();
  auto diag_f = diag->get_diagnostic();

  // Check result: diag should be filled with fill_value
  auto some_field = qc.clone();
  some_field.deep_copy(constants::fill_value<Real>);
  REQUIRE(views_are_equal(diag_f, some_field));

  const Real a_day = 24.0 * 60.0 * 60.0;  // seconds

  constexpr int ntests = 10;

  for(int itest = 2; itest < ntests; itest++) {
    // Run diag again
    some_field.deep_copy(qc);

    diag->init_timestep(t0);

    util::TimeStamp t1({2024, 1, itest}, {0, 0, 0});  // a day later
    qc.get_header().get_tracking().update_time_stamp(t1);
    randomize_uniform(qc, seed++, 0, 200);

    diag->compute_diagnostic();
    some_field.update(qc, 1.0 / a_day, -1.0 / a_day);
    REQUIRE(views_are_equal(diag_f, some_field));

    // reset t0 to t1 to keep iterating...
    t0 = t1;
  }
}

}  // namespace scream
