#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

TEST_CASE("field_prev") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp ts({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager
  const int nlevs = 7;
  const int ncols = 2;

  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);

  // Input (randomized) qc
  FieldLayout layout = grid->get_3d_scalar_layout(LEV);
  FieldIdentifier qc_fid("qc", layout, kg / kg, grid->name());

  Field qc(qc_fid);
  qc.allocate_view();

  // Random number generator seed
  int seed = get_random_test_seed(&comm);

  // Construct the Diagnostics
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  // Create and set up the diagnostic
  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("FieldPrev", comm, params, grid)); // Bad construction

  params.set("grid_name", grid->name());
  REQUIRE_THROWS(diag_factory.create("FieldPrev", comm, params, grid)); // Still no field_name

  // Set time for qc and randomize its values
  randomize_uniform(qc, seed++, 0, 200);

  // Create and set up the diagnostic
  params.set<std::string>("field_name", "qc");
  auto diag = diag_factory.create("FieldPrev", comm, params, grid);
  diag->set_input_field(qc);
  diag->initialize();

  REQUIRE (diag->get_diagnostic().has_valid_mask());

  // Simulate a derived diagnostic whose source has no valid timestamp at
  // init_timestep time (i.e. it hasn't been computed yet on step 1).
  diag->init_timestep(ts);

  // Source is now "computed" (gets a valid timestamp and random values), but it's too late...
  qc.get_header().get_tracking().update_time_stamp(ts);
  randomize_uniform(qc, seed++, 0, 200);
  diag->compute_diagnostic(ts);
  auto diag_mask = diag->get_diagnostic().get_valid_mask();

  auto tgt_mask = diag_mask.clone();
  tgt_mask.deep_copy(0);

  REQUIRE(views_are_equal(diag_mask,tgt_mask));

  constexpr int ntests = 10;
  tgt_mask.deep_copy(1);
  for (int itest = 2; itest < ntests; itest++) {
    // Save current qc before we advance
    auto qc_prev = qc.clone();
    qc_prev.deep_copy(qc);

    // init_timestep stores the current qc as the "previous" value
    diag->init_timestep(ts);

    // Advance time and update qc
    ts += 300;
    qc.get_header().get_tracking().update_time_stamp(ts);
    randomize_uniform(qc, seed++, 0, 200);

    // Diagnostic should return the value stored at init_timestep (qc_prev)
    diag->compute_diagnostic(ts);
    auto diag_f = diag->get_diagnostic();
    REQUIRE(views_are_equal(diag_f, qc_prev));

    REQUIRE(views_are_equal(diag_mask,tgt_mask));
  }
}

}  // namespace scream
