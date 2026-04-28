#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {


TEST_CASE("field_over_dt") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);

  util::TimeStamp ts({2024, 1, 1}, {0, 0, 0});

  const int nlevs = 7;
  const int ncols = 2;

  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);

  FieldLayout scalar2d_layout = grid->get_3d_scalar_layout(LEV);
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());

  Field qc(qc_fid);
  qc.allocate_view();

  int seed = get_random_test_seed(&comm);

  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("FieldOverDt", comm,
                                     params, grid));  // No 'field_name'

  // Set time for qc and randomize its values
  qc.get_header().get_tracking().update_time_stamp(ts);
  randomize_uniform(qc, seed++, 0, 200);

  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("field_name", "qc");
  auto diag = diag_factory.create("FieldOverDt", comm, params, grid);
  diag->set_input_field(qc);
  diag->initialize();

  constexpr int ntests = 10;
  const Real a_day = 24.0 * 60.0 * 60.0;  // seconds

  auto diag_f = diag->get_diagnostic();
  for (int itest = 2; itest < ntests; itest++) {
    // Save qc before advancing
    auto qc_prev = qc.clone();
    qc_prev.deep_copy(qc);

    // init_timestep saves start of step
    diag->init_timestep(ts);

    // Advance time by one day and update qc
    ts += a_day;
    qc.get_header().get_tracking().update_time_stamp(ts);
    randomize_uniform(qc, seed++, 0, 200);
    qc.sync_to_dev();

    // Diagnostic should be qc(ts) / dt
    diag->compute_diagnostic(ts);
    diag_f.sync_to_host();

    auto qc_v    = qc.get_view<Real**, Host>();
    auto diag_v  = diag_f.get_view<Real**, Host>();
    for (int icol = 0; icol < ncols; ++icol)
      for (int ilev = 0; ilev < nlevs; ++ilev)
        REQUIRE(diag_v(icol, ilev) == Approx(qc_v(icol, ilev) / a_day));
  }
}

}  // namespace scream
