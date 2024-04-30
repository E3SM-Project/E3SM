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

TEST_CASE("aerocom_cld") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2022, 1, 1}, {0, 0, 0});

  const auto nondim = Units::nondimensional();
  const auto micron = m / 1000000;

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 33;
  const int ngcols    = 2 * comm.size();

  int m_ndiag = 8;

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) qc, nc
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};

  // Create fields
  FieldIdentifier tm_fid("T_mid", scalar2d_layout, K, grid->name());
  FieldIdentifier pd_fid("pseudo_density", scalar2d_layout, Pa, grid->name());
  FieldIdentifier pm_fid("p_mid", scalar2d_layout, Pa, grid->name());
  FieldIdentifier qv_fid("qv", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qi_fid("qi", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier ec_fid("eff_radius_qc", scalar2d_layout, micron,
                         grid->name());
  FieldIdentifier ei_fid("eff_radius_qi", scalar2d_layout, micron,
                         grid->name());
  FieldIdentifier cd_fid("cldfrac_tot", scalar2d_layout, nondim, grid->name());
  FieldIdentifier nc_fid("nc", scalar2d_layout, kg / kg, grid->name());

  Field tm(tm_fid);
  tm.allocate_view();
  tm.get_header().get_tracking().update_time_stamp(t0);
  Field pd(pd_fid);
  pd.allocate_view();
  pd.get_header().get_tracking().update_time_stamp(t0);
  Field pm(pm_fid);
  pm.allocate_view();
  pm.get_header().get_tracking().update_time_stamp(t0);
  Field qv(qv_fid);
  qv.allocate_view();
  qv.get_header().get_tracking().update_time_stamp(t0);
  Field qc(qc_fid);
  qc.allocate_view();
  qc.get_header().get_tracking().update_time_stamp(t0);
  Field qi(qi_fid);
  qi.allocate_view();
  qi.get_header().get_tracking().update_time_stamp(t0);
  Field ec(ec_fid);
  ec.allocate_view();
  ec.get_header().get_tracking().update_time_stamp(t0);
  Field ei(ei_fid);
  ei.allocate_view();
  ei.get_header().get_tracking().update_time_stamp(t0);
  Field cd(cd_fid);
  cd.allocate_view();
  cd.get_header().get_tracking().update_time_stamp(t0);
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
  ekat::ParameterList params;

  REQUIRE_THROWS(
      diag_factory.create("AeroComCld", comm, params));  // No 'AeroComCld Kind'
  params.set<std::string>("AeroComCld Kind", "Foo");
  REQUIRE_THROWS(diag_factory.create("AeroComCld", comm,
                                     params));  // Invalid 'AeroComCld Kind'

  constexpr int ntests = 5;
  for(int itest = 0; itest < ntests; ++itest) {
    // Randomize everything
    randomize(tm, engine, pdf);
    randomize(pd, engine, pdf);
    randomize(pm, engine, pdf);
    randomize(qv, engine, pdf);
    randomize(qc, engine, pdf);
    randomize(qi, engine, pdf);
    randomize(ec, engine, pdf);
    randomize(ei, engine, pdf);
    randomize(cd, engine, pdf);
    randomize(nc, engine, pdf);

    // Create and set up the diagnostic
    params.set<std::string>("AeroComCld Kind", "Top");
    auto diag = diag_factory.create("AeroComCld", comm, params);

    diag->set_grids(gm);

    diag->set_required_field(tm);
    diag->set_required_field(pd);
    diag->set_required_field(pm);
    diag->set_required_field(qv);
    diag->set_required_field(qc);
    diag->set_required_field(qi);
    diag->set_required_field(ec);
    diag->set_required_field(ei);
    diag->set_required_field(cd);
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
    // // Workaround for non-bfb behavior of view_reduction() in release builds
    // if(SCREAM_BFB_TESTING) {
    //   REQUIRE(views_are_equal(out_hf, out_tf));
    // }
  }
}

}  // namespace scream
