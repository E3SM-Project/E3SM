#include "catch2/catch.hpp"
#include "diagnostics/aerocom_cld_util.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

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
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  const auto nondim = Units::nondimensional();
  const auto micron = m / 1000000;

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 9;
  const int ngcols    = 1 * comm.size();

  // See how many diags we are calculating
  AeroComCldDiagUtil aercom_util;
  int ndiags = aercom_util.size;

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input
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
  FieldIdentifier nc_fid("nc", scalar2d_layout, 1 / kg, grid->name());
  FieldIdentifier ni_fid("ni", scalar2d_layout, 1 / kg, grid->name());

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
  Field ni(ni_fid);
  ni.allocate_view();
  ni.get_header().get_tracking().update_time_stamp(t0);

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

  constexpr int ntests = 3;
  for(int itest = 0; itest < ntests; ++itest) {
    // Randomize everything to add ensure resiliency
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
    randomize(ni, engine, pdf);

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
    diag->set_required_field(ni);

    diag->initialize(t0, RunType::Initial);

    // Case 1: if the cloud fraction is zero, everything is zero
    cd.deep_copy(0.0);
    diag->compute_diagnostic();
    Field diag_f = diag->get_diagnostic();
    diag_f.sync_to_host();
    auto diag_v = diag_f.get_view<Real **, Host>();
    for(int idiag = 0; idiag < ndiags; ++idiag) {
      REQUIRE(diag_v(0, idiag) == Real(0.0));
    }

    // Case 2: if the cloud fraction is one, everything takes 1 * its value
    // Take a moment to set "sensible" for other things...
    cd.deep_copy(1.0);
    tm.deep_copy(300.0);
    pd.deep_copy(10.0);
    pm.deep_copy(100.0);
    qv.deep_copy(1.0);
    qc.deep_copy(1.0);
    qi.deep_copy(1.0);
    ec.deep_copy(10.0);
    ei.deep_copy(10.0);
    nc.deep_copy(5.0);
    ni.deep_copy(1.0);
    diag->compute_diagnostic();
    diag->get_diagnostic().sync_to_host();
    diag_f = diag->get_diagnostic();
    diag_v = diag_f.get_view<Real **, Host>();
    for(int icol = 0; icol < grid->get_num_local_dofs(); ++icol) {
      REQUIRE(diag_v(icol, 0) == Real(300.0));
      REQUIRE(diag_v(icol, 1) == Real(100.0));
      REQUIRE(diag_v(icol, 2) == Real(0.5));
      REQUIRE(diag_v(icol, 3) == Real(0.5));
      REQUIRE(diag_v(icol, 4) > Real(0.0));
      REQUIRE(diag_v(icol, 5) > Real(0.0));
      REQUIRE(diag_v(icol, 6) > Real(0.0));
      REQUIRE(diag_v(icol, 7) == Real(1.0));
    }

    // Case 3: test the max overlap (if contiguous cloudy layers, then max)
    cd.deep_copy<Host>(0);
    auto cd_v  = cd.get_view<Real **, Host>();
    cd_v(0, 1) = 0.5;
    cd_v(0, 2) = 0.7;  // ------> max!
    cd_v(0, 3) = 0.3;
    cd_v(0, 4) = 0.2;
    cd.sync_to_dev();
    diag->compute_diagnostic();
    diag_f = diag->get_diagnostic();
    diag_f.sync_to_host();
    diag_v = diag_f.get_view<Real **, Host>();
    REQUIRE(diag_v(0, 7) == Real(0.7));

    // Case 3xtra: test max overlap again
    // This case should produce >0.7 due to slight enhancement in the presence
    // of a local minimum (0.1 is the local minimum between 0.2 and 0.4)
    cd_v(0, 5) = 0.1;
    cd_v(0, 6) = 0.4;
    cd_v(0, 7) = 0.2;
    cd.sync_to_dev();
    diag->compute_diagnostic();
    diag_f = diag->get_diagnostic();
    diag_f.sync_to_host();
    diag_v = diag_f.get_view<Real **, Host>();
    REQUIRE(diag_v(0, 7) > Real(0.7));

    // Case 4: test random overlap
    // If non-contiguous cloudy layers, then random
    cd_v(0, 6) = 0.0;
    cd_v(0, 7) = 0.1;
    cd.sync_to_dev();
    diag->compute_diagnostic();
    diag_f = diag->get_diagnostic();
    diag_f.sync_to_host();
    diag_v = diag_f.get_view<Real **, Host>();
    REQUIRE(diag_v(0, 7) > Real(0.7));  // must be larger than the max!

    // Case 5a: test independence of ice and liq fractions
    cd_v(0, 3) = 1.0;
    cd_v(0, 7) = 1.0;
    cd_v(0, 8) = 0.2;
    qc.deep_copy(1.0);
    qi.deep_copy(0.0);  // zero ice!
    cd.sync_to_dev();
    diag->compute_diagnostic();
    diag_f = diag->get_diagnostic();
    diag_f.sync_to_host();
    diag_v = diag_f.get_view<Real **, Host>();
    REQUIRE(diag_v(0, 7) == Real(1.0));
    REQUIRE(diag_v(0, 3) == Real(1.0));
    REQUIRE(diag_v(0, 2) == Real(0.0));  // zero ice!

    // Case 5b: test independence of ice and liq fractions
    qc.deep_copy(0.0);  // zero liq!
    qi.deep_copy(1.0);
    diag->compute_diagnostic();
    diag_f = diag->get_diagnostic();
    diag_f.sync_to_host();
    diag_v = diag_f.get_view<Real **, Host>();
    REQUIRE(diag_v(0, 7) == Real(1.0));
    REQUIRE(diag_v(0, 3) == Real(0.0));  // zero liq!
    REQUIRE(diag_v(0, 2) == Real(1.0));

    // Case 6: test independence of ice and liquid fractions
    // There is NOT complete independence...
    // Essentially, higher ice clouds mask lower liquid clouds
    // This can be problematic if the ice clouds are thin...
    // We will revisit and validate this assumption later
    auto qc_v = qc.get_view<Real **, Host>();
    auto qi_v = qi.get_view<Real **, Host>();
    cd.deep_copy<Host>(0);
    cd_v(0, 1) = 0.5;  // ice
    cd_v(0, 2) = 0.7;  // ice ------> max!
    cd_v(0, 3) = 0.3;  // ice
    // note cd_v(0, 4) is 0.0
    cd_v(0, 5) = 0.2;  // liq
    cd_v(0, 6) = 0.5;  // liq ------> not max!
    cd_v(0, 7) = 0.1;  // liq
    // note cd_v(0, 8) is 0.0
    qi.deep_copy<Host>(0);
    qi_v(0, 1) = 100;
    qi_v(0, 2) = 200;
    qi_v(0, 3) = 50;
    // note qi_v(0, 4) = 0.0
    qc.deep_copy<Host>(0);
    // note qc_v(0, 4) = 0.0
    qc_v(0, 5) = 20;
    qc_v(0, 6) = 50;
    qc_v(0, 7) = 10;
    cd.sync_to_dev();
    qi.sync_to_dev();
    qc.sync_to_dev();
    diag->compute_diagnostic();
    diag_f = diag->get_diagnostic();
    diag_f.sync_to_host();
    diag_v = diag_f.get_view<Real **, Host>();
    REQUIRE(diag_v(0, 7) > Real(0.7));   // unaffected (see test case 4)
    REQUIRE(diag_v(0, 3) < Real(0.5));   // not max!
    REQUIRE(diag_v(0, 2) == Real(0.7));  // max!
  }
}

}  // namespace scream
