#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

TEST_CASE("aodvis") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2022, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  const int nlevs = 33;
  const int ncols = 10;

  int nbnds = eamxx_swbands();
  int swvis = eamxx_vis_swband_idx();

  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);

  // Input (randomized) tau
  FieldLayout scalar3d_swband_layout =
      grid->get_3d_vector_layout(LEV, nbnds, "swband");
  FieldIdentifier tau_fid("aero_tau_sw", scalar3d_swband_layout, none, grid->name());
  Field tau(tau_fid);
  tau.allocate_view();
  tau.get_header().get_tracking().update_time_stamp(t0);
  // Input (randomized) sunlit
  FieldLayout scalar2d_layout = grid->get_2d_scalar_layout();
  FieldIdentifier sunlit_fid("sunlit_mask", scalar2d_layout, none, grid->name(), DataType::IntType);
  Field sunlit(sunlit_fid);
  sunlit.allocate_view();
  sunlit.get_header().get_tracking().update_time_stamp(t0);

  // Random number generator seed
  int seed = get_random_test_seed(&comm);

  // Construct the Diagnostics
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  constexpr int ntests = 5;
  for(int itest = 0; itest < ntests; ++itest) {
    // Randomize tau
    randomize_uniform(tau, seed++, 0, 0.005);

    // Randomize sunlit
    randomize_uniform(sunlit, seed++, 0, 1);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("AerosolOpticalDepth550nm", comm, params, grid);
    diag->set_input_field(tau);
    diag->set_input_field(sunlit);
    diag->initialize();

    // Run diag
    diag->compute_diagnostic(t0);

    // Check result
    tau.sync_to_host();
    diag->get_diagnostic().sync_to_host();

    const auto sun_h  = sunlit.get_view<const int *, Host>();
    const auto tau_h  = tau.get_view<const Real ***, Host>();
    const auto aod_hf = diag->get_diagnostic();
    const auto aod_mask = aod_hf.get_valid_mask();

    Field aod_tf = diag->get_diagnostic().clone();
    auto aod_t = aod_tf.get_view<Real *, Host>();

    for(int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
      if(sun_h(icol)) {
        aod_t(icol) = 0;
        for(int ilev = 0; ilev < nlevs; ++ilev) {
          aod_t(icol) += tau_h(icol, swvis, ilev);
        }
      } else {
        aod_t(icol) = constants::fill_value<Real>;
      }
    }
    aod_hf.sync_to_dev();
    aod_tf.sync_to_dev();

    // Workaround for non-bfb behavior of view_reduction() in release builds
    if(SCREAM_BFB_TESTING) {
      REQUIRE(views_are_equal(aod_hf, aod_tf));
    }
    // Ensure the masks are identical (should be pointing exactly to each other)
    REQUIRE(views_are_equal(aod_mask, sunlit));
  }
}

}  // namespace scream
