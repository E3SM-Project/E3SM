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

TEST_CASE("aerocom_cldtop") {
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

  int m_ndiag = 2;

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) qc, nc
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier nc_fid("nc", scalar2d_layout, kg / kg, grid->name());
  Field qc(qc_fid);
  qc.allocate_view();
  qc.get_header().get_tracking().update_time_stamp(t0);
  Field nc(nc_fid);
  nc.allocate_view();
  nc.get_header().get_tracking().update_time_stamp(t0);

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
    // Randomize qc, nc
    randomize(qc, engine, pdf);
    randomize(nc, engine, pdf);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("AeroComCldTop", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(qc);
    diag->set_required_field(nc);
    diag->initialize(t0, RunType::Initial);

    // Run diag
    diag->compute_diagnostic();

    // Check result
    qc.sync_to_host();
    nc.sync_to_host();
    diag->get_diagnostic().sync_to_host();

    const auto qc_h   = qc.get_view<const Real **, Host>();
    const auto nc_h   = nc.get_view<const Real **, Host>();
    const auto out_hf = diag->get_diagnostic();

    Field out_tf = diag->get_diagnostic().clone();
    out_tf.deep_copy<double, Host>(0.0);
    auto out_t = out_tf.get_view<Real **, Host>();

    for(int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
      for(int ilev = 0; ilev < nlevs; ++ilev) {
        out_t(icol, 0) += qc_h(icol, ilev);
        out_t(icol, 1) += nc_h(icol, ilev);
      }
    }
    out_hf.sync_to_dev();
    out_tf.sync_to_dev();
    // Workaround for non-bfb behavior of view_reduction() in release builds
    if(SCREAM_BFB_TESTING) {
      REQUIRE(views_are_equal(out_hf, out_tf));
    }
  }
}

}  // namespace scream
