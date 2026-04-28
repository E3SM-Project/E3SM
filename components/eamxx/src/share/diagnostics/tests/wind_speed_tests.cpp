#include "catch2/catch.hpp"

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

namespace scream {


TEST_CASE("wind_speed")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager - single column for these tests
  const int nlevs = 33;
  const int ncols = 2;;
  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);

  // Input (randomized) velocity
  auto vector3d = grid->get_3d_vector_layout(LEV,2,"dim");
  FieldIdentifier uv_fid ("horiz_winds",vector3d,m/s,grid->name());
  Field uv(uv_fid);
  uv.allocate_view();
  uv.get_header().get_tracking().update_time_stamp(t0);

  // Random number generator seed
  int seed = get_random_test_seed(&comm);

  // Construct the Diagnostics
  auto& diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  constexpr int ntests = 5;
#ifdef NDEBUG
  constexpr int ulp_tol = 1;
#else
  constexpr int ulp_tol = 0;
#endif
  for (int itest=0; itest<ntests; ++itest) {
    // Randomize wind
    randomize_uniform(uv,seed++);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("wind_speed",comm,params, grid);
    diag->set_input_field(uv);
    diag->initialize();

    // Run diag
    diag->compute_diagnostic(t0);

    // Check result
    uv.sync_to_host();
    diag->get_diagnostic().sync_to_host();

    auto uv_h = uv.get_view<const Real***,Host>();
    auto ws_h = diag->get_diagnostic().get_view<const Real**,Host>();

    for (int icol=0; icol<grid->get_num_local_dofs(); ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        const auto u = uv_h (icol,0,ilev);
        const auto v = uv_h (icol,1,ilev);
        REQUIRE_THAT (ws_h(icol,ilev), Catch::Matchers::WithinULP(std::sqrt(u*u+v*v),ulp_tol));
      }
    }
  }
}

} // namespace
