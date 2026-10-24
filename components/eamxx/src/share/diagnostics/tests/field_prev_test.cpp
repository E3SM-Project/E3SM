#include <catch2/catch.hpp>

#include "share/util/eamxx_universal_constants.hpp"
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

  REQUIRE (diag->get().has_valid_mask());

  // Simulate a derived diagnostic whose source has no valid timestamp at
  // init_timestep time (i.e. it hasn't been computed yet on step 1).
  diag->init_timestep(ts);

  // Source is now "computed" (gets a valid timestamp and random values), but it's too late...
  qc.get_header().get_tracking().update_time_stamp(ts);
  randomize_uniform(qc, seed++, 0, 200);
  diag->compute(ts);
  auto diag_mask = diag->get().get_valid_mask();

  auto tgt_mask = diag_mask.clone();
  tgt_mask.deep_copy(0);

  REQUIRE(views_are_equal(diag_mask,tgt_mask));

  constexpr int ntests = 10;
  tgt_mask.deep_copy(1);
  for (int itest = 2; itest < ntests; itest++) {
    // Save current qc before we advance
    auto qc_prev = qc.clone(CloneFlags::CopyData);

    // init_timestep stores the current qc as the "previous" value
    diag->init_timestep(ts);

    // Advance time and update qc
    ts += 300;
    qc.get_header().get_tracking().update_time_stamp(ts);
    randomize_uniform(qc, seed++, 0, 200);

    // Diagnostic should return the value stored at init_timestep (qc_prev)
    diag->compute(ts);
    auto diag_f = diag->get();
    REQUIRE(views_are_equal(diag_f, qc_prev));

    REQUIRE(views_are_equal(diag_mask,tgt_mask));
  }
}

TEST_CASE("field_prev_masked") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  util::TimeStamp ts({2024, 1, 1}, {0, 0, 0});

  const int nlevs = 7;
  const int ncols = 2;
  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);
  const int nlocal = grid->get_num_local_dofs();

  FieldIdentifier qc_fid("qc", grid->get_3d_scalar_layout(LEV), kg / kg, grid->name());
  Field qc(qc_fid,true);
  qc.create_valid_mask(Field::MaskInit::Valid);

  int seed = get_random_test_seed(&comm);

  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  params.set("grid_name", grid->name());
  params.set<std::string>("field_name", "qc");
  auto diag = diag_factory.create("FieldPrev", comm, params, grid);
  diag->set_input_field(qc);
  diag->initialize();

  constexpr auto fv = constants::fill_value<Real>;

  // Step 1: qc is valid everywhere
  qc.get_header().get_tracking().update_time_stamp(ts);
  randomize_uniform(qc, seed++, 1, 200);
  diag->init_timestep(ts);
  ts += 300;
  diag->compute(ts);

  auto d = diag->get();
  auto d_mask = d.get_valid_mask();
  d.sync_to_host();
  d_mask.sync_to_host();
  auto d_v = d.get_view<const Real**,Host>();
  auto m_v = d_mask.get_view<const int**,Host>();
  for (int i=0; i<nlocal; ++i) {
    for (int k=0; k<nlevs; ++k) {
      REQUIRE (m_v(i,k)==1);
      REQUIRE (d_v(i,k)!=fv);
    }
  }

  // Step 2: level 0 of qc turns invalid, so the diag must hold fill_value there
  {
    auto qm = qc.get_valid_mask();
    qm.sync_to_host();
    auto qm_h = qm.get_view<int**,Host>();
    for (int i=0; i<nlocal; ++i) qm_h(i,0) = 0;
    qm.sync_to_dev();
  }
  randomize_uniform(qc, seed++, 1, 200);
  qc.get_header().get_tracking().update_time_stamp(ts);
  auto qc_prev = qc.clone(CloneFlags::CopyData);
  qc_prev.sync_to_host();
  auto p_v = qc_prev.get_view<const Real**,Host>();
  diag->init_timestep(ts);
  ts += 300;
  diag->compute(ts);

  d.sync_to_host();
  d_mask.sync_to_host();
  for (int i=0; i<nlocal; ++i) {
    for (int k=0; k<nlevs; ++k) {
      if (k==0) {
        REQUIRE (m_v(i,k)==0);
        REQUIRE (d_v(i,k)==fv);
      } else {
        REQUIRE (m_v(i,k)==1);
        REQUIRE (d_v(i,k)==p_v(i,k));
      }
    }
  }
}

}  // namespace scream
