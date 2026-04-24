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

TEST_CASE("field_prev") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager
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
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  // Create and set up the diagnostic
  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("FieldPrevDiag", comm, params)); // Bad construction

  params.set("grid_name", grid->name());
  REQUIRE_THROWS(diag_factory.create("FieldPrevDiag", comm, params)); // Still no field_name

  // Set time for qc and randomize its values
  qc.get_header().get_tracking().update_time_stamp(t0);
  randomize_uniform(qc, seed++, 0, 200);

  // Create and set up the diagnostic
  params.set<std::string>("field_name", "qc");
  auto diag = diag_factory.create("FieldPrevDiag", comm, params);
  diag->set_grid(grid);
  diag->set_required_field(qc);
  diag->initialize();

  // Run diag before any init_timestep call: qc already has a valid t0
  // timestamp, so the fallback path returns X(t=0) = qc rather than fill_value.
  diag->compute_diagnostic();
  auto diag_f = diag->get_diagnostic();
  REQUIRE(views_are_equal(diag_f, qc));

  // Simulate a derived diagnostic whose source has no valid timestamp at
  // init_timestep time (i.e. it hasn't been computed yet on step 1).
  // After init_timestep the source is computed and gets a valid timestamp;
  // the diagnostic should then return the current source value (X_prev = X(t=0))
  // rather than fill_value or stale zeros — this is the key test for the
  // "first step of a derived _prev field" correctness.
  {
    FieldIdentifier qc_uninit_fid("qc_uninit", scalar2d_layout, kg / kg, grid->name());
    Field qc_uninit(qc_uninit_fid);
    qc_uninit.allocate_view();  // zero-initialized, no valid timestamp

    ekat::ParameterList p2;
    p2.set("grid_name", grid->name());
    p2.set<std::string>("field_name", "qc_uninit");
    auto diag2 = diag_factory.create("FieldPrevDiag", comm, p2);
    diag2->set_grid(grid);
    diag2->set_required_field(qc_uninit);
    diag2->initialize();

    // init_timestep: source has no valid timestamp → no capture
    diag2->init_timestep(t0);

    // Source is now "computed" (gets a valid timestamp and random values),
    // simulating what the output manager does before calling this diagnostic.
    qc_uninit.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(qc_uninit, seed++, 0, 200);

    // Fallback: m_f_prev invalid, source now valid → returns current source value
    diag2->compute_diagnostic();
    auto diag2_f = diag2->get_diagnostic();
    REQUIRE(views_are_equal(diag2_f, qc_uninit));
  }

  constexpr int ntests = 10;

  for (int itest = 2; itest < ntests; itest++) {
    // Save current qc before we advance
    auto qc_prev = qc.clone();
    qc_prev.deep_copy(qc);

    // init_timestep stores the current qc as the "previous" value
    diag->init_timestep(t0);

    // Advance time and update qc
    util::TimeStamp t1({2024, 1, itest}, {0, 0, 0});
    qc.get_header().get_tracking().update_time_stamp(t1);
    randomize_uniform(qc, seed++, 0, 200);

    // Diagnostic should return the value stored at init_timestep (qc_prev)
    diag->compute_diagnostic();
    REQUIRE(views_are_equal(diag_f, qc_prev));

    t0 = t1;
  }
}

}  // namespace scream
