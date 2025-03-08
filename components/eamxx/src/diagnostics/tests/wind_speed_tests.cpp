#include "catch2/catch.hpp"

#include "diagnostics/register_diagnostics.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

std::shared_ptr<GridsManager>
create_gm (const ekat::Comm& comm, const int ncols, const int nlevs) {

  const int num_global_cols = ncols*comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names",vos_t{"Point Grid"});
  auto& pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type","point_grid");
  pl.set("aliases",vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm,gm_params);
  gm->build_grids();

  return gm;
}

TEST_CASE("wind_speed")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0 ({2022,1,1},{0,0,0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 33;
  const int ngcols = 2*comm.size();;
  auto gm = create_gm(comm,ngcols,nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) velocity
  auto vector3d = grid->get_3d_vector_layout(true,2);
  FieldIdentifier uv_fid ("horiz_winds",vector3d,m/s,grid->name());
  Field uv(uv_fid);
  uv.allocate_view();
  uv.get_header().get_tracking().update_time_stamp(t0);

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(-1,1);
  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string,std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  constexpr int ntests = 5;
#ifdef NDEBUG
  constexpr int ulp_tol = 1;
#else
  constexpr int ulp_tol = 0;
#endif
  for (int itest=0; itest<ntests; ++itest) {
    // Randomize wind
    randomize(uv,engine,pdf);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("wind_speed",comm,params);
    diag->set_grids(gm);
    diag->set_required_field(uv);
    diag->initialize(t0,RunType::Initial);

    // Run diag
    diag->compute_diagnostic();

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
