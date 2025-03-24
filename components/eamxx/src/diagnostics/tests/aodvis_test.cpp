#include <share/field/field_utils.hpp>

#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"
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

  Real var_fill_value = constants::DefaultFillValue<Real>().value;

  Real some_limit = 0.0025;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2022, 1, 1}, {0, 0, 0});

  const auto nondim = Units::nondimensional();

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 33;
  const int ngcols    = 10 * comm.size();

  int nbnds = eamxx_swbands();
  int swvis = eamxx_vis_swband_idx();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) tau
  FieldLayout scalar3d_swband_layout =
      grid->get_3d_vector_layout(true, nbnds, "swband");
  FieldIdentifier tau_fid("aero_tau_sw", scalar3d_swband_layout, nondim,
                          grid->name());
  Field tau(tau_fid);
  tau.allocate_view();
  tau.get_header().get_tracking().update_time_stamp(t0);
  // Input (randomized) sunlit
  FieldLayout scalar2d_layout = grid->get_2d_scalar_layout();
  FieldIdentifier sunlit_fid("sunlit", scalar2d_layout, nondim, grid->name());
  Field sunlit(sunlit_fid);
  sunlit.allocate_view();
  sunlit.get_header().get_tracking().update_time_stamp(t0);

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0, 0.005);
  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  constexpr int ntests = 5;
  for(int itest = 0; itest < ntests; ++itest) {
    // Randomize tau
    randomize(tau, engine, pdf);

    // Randomize sunlit
    randomize(sunlit, engine, pdf);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("AerosolOpticalDepth550nm", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(tau);
    diag->set_required_field(sunlit);
    diag->initialize(t0, RunType::Initial);

    auto sun_h = sunlit.get_view<Real *, Host>();
    for(int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
      // zero out all sun_h if below 0.05
      if(sun_h(icol) < some_limit) {
        sun_h(icol) = 0.0;
      }
    }
    // sync to device
    sunlit.sync_to_dev();

    // Run diag
    diag->compute_diagnostic();

    // Check result
    tau.sync_to_host();
    diag->get_diagnostic().sync_to_host();

    const auto tau_h  = tau.get_view<const Real ***, Host>();
    const auto aod_hf = diag->get_diagnostic();

    Field aod_tf = diag->get_diagnostic().clone();
    aod_tf.deep_copy<double, Host>(0.0);
    auto aod_t = aod_tf.get_view<Real *, Host>();

    for(int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
      if(sun_h(icol) < some_limit) {
        aod_t(icol) = var_fill_value;
      } else {
        for(int ilev = 0; ilev < nlevs; ++ilev) {
          aod_t(icol) += tau_h(icol, swvis, ilev);
        }
      }
    }
    aod_hf.sync_to_dev();
    aod_tf.sync_to_dev();

    // Workaround for non-bfb behavior of view_reduction() in release builds
    if(SCREAM_BFB_TESTING) {
      REQUIRE(views_are_equal(aod_hf, aod_tf));
    }
  }
}

}  // namespace scream
