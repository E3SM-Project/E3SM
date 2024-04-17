#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/scream_setup_random_test.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols,
                                        const int nlevs) {
  const int num_global_cols = ncols * comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names", vos_t{"Point Grid"});
  auto &pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type", "point_grid");
  pl.set("aliases", vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", num_global_cols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm, gm_params);
  gm->build_grids();

  return gm;
}

TEST_CASE("aodvis") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2022, 1, 1}, {0, 0, 0});

  const auto nondim = Units::nondimensional();

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 33;
  const int ngcols    = 2 * comm.size();

  // TODO: Don't hardcode this!
  constexpr int nbnds = 14;

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) tau
  FieldLayout scalar3d_swband_layout{{COL, SWBND, LEV}, {ngcols, nbnds, nlevs}};
  FieldIdentifier tau_fid("aero_tau_sw", scalar3d_swband_layout, nondim,
                          grid->name());
  Field tau(tau_fid);
  tau.allocate_view();
  tau.get_header().get_tracking().update_time_stamp(t0);

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0, 0.05);
  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  constexpr int ntests = 5;
  for(int itest = 0; itest < ntests; ++itest) {
    // Randomize tau
    randomize(tau, engine, pdf);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("AerosolOpticalDepth550nm", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(tau);
    diag->initialize(t0, RunType::Initial);

    // Run diag
    diag->compute_diagnostic();

    // Check result
    tau.sync_to_host();
    diag->get_diagnostic().sync_to_host();

    auto tau_h  = tau.get_view<const Real ***, Host>();
    auto aod_hf = diag->get_diagnostic();
    auto aod_h  = aod_hf.get_view<const Real *, Host>();

    Field aod_tf = diag->get_diagnostic().clone();
    aod_tf.deep_copy<double, Host>(0.0);
    auto aod_t = aod_tf.get_view<Real *, Host>();

    for(int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
      auto aod_temp = 0.0;
      for(int ilev = 0; ilev < nlevs; ++ilev) {
        // TODO: Don't hardcode the 10!
        aod_temp += tau_h(icol, 10, ilev);
        aod_t(icol) += tau_h(icol, 10, ilev);
      }
      // Allow for sp tolerance?
      REQUIRE(std::abs(aod_temp - aod_h(icol)) < 1e-6);
    }
    // Another way to test is to compare the views
    REQUIRE(views_are_equal(aod_hf, aod_tf));
  }
}

}  // namespace scream
