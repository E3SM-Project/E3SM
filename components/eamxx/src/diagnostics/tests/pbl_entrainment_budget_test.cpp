#include "catch2/catch.hpp"
#include "diagnostics/pbl_entrainment_budget_util.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include <ekat_test_utils.hpp>
#include "share/field/field_utils.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
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

TEST_CASE("pbl_entrainment_budget") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using PF = scream::PhysicsFunctions<DefaultDevice>;
  using PC = scream::physics::Constants<Real>;

  // A numerical tolerance
  // Accumulations in the Kokkos threaded reductions may be done in a
  // different order than the manual ones below, so we can only test
  // correctness up to a tolerance
  auto tol = std::numeric_limits<Real>::epsilon() * 100;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  scorpio::init_subsystem(comm);

  std::string iop_file = ekat::TestSession::get().params.at("iop_file");

  scorpio::register_file(iop_file, scorpio::FileMode::Read);
  const int nlevs = scorpio::get_dimlen(iop_file, "lev");
  // From now on, when we read vars, "time" must be treated as unlimited, to
  // avoid issues
  if(not scorpio::has_time_dim(iop_file)) {
    scorpio::mark_dim_as_time(iop_file, "time");
  }

  // time stamps
  util::TimeStamp t0({2022, 1, 1}, {0, 1, 0});
  util::TimeStamp t1({2022, 1, 1}, {0, 2, 0});
  util::TimeStamp t2({2022, 1, 1}, {0, 3, 0});

  // Create a grids manager - single column for these tests
  const int ngcols = 1 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qv_fid("qv", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier tm_fid("T_mid", scalar2d_layout, K, grid->name());
  FieldIdentifier pm_fid("p_mid", scalar2d_layout, Pa, grid->name());
  FieldIdentifier pd_fid("pseudo_density", scalar2d_layout, Pa, grid->name());
  FieldIdentifier sd_fid("SW_flux_dn", scalar2d_layout, W / m * m,
                         grid->name());
  FieldIdentifier su_fid("SW_flux_up", scalar2d_layout, W / m * m,
                         grid->name());
  FieldIdentifier ld_fid("LW_flux_dn", scalar2d_layout, W / m * m,
                         grid->name());
  FieldIdentifier lu_fid("LW_flux_up", scalar2d_layout, W / m * m,
                         grid->name());

  Field qc(qc_fid);
  qc.allocate_view();
  qc.get_header().get_tracking().update_time_stamp(t0);
  Field qv(qv_fid);
  qv.allocate_view();
  qv.get_header().get_tracking().update_time_stamp(t0);
  Field tm(tm_fid);
  tm.allocate_view();
  tm.get_header().get_tracking().update_time_stamp(t0);
  Field pm(pm_fid);
  pm.allocate_view();
  pm.get_header().get_tracking().update_time_stamp(t0);
  Field pd(pd_fid);
  pd.allocate_view();
  pd.get_header().get_tracking().update_time_stamp(t0);
  Field sd(sd_fid);
  sd.allocate_view();
  sd.get_header().get_tracking().update_time_stamp(t0);
  Field su(su_fid);
  su.allocate_view();
  su.get_header().get_tracking().update_time_stamp(t0);
  Field ld(ld_fid);
  ld.allocate_view();
  ld.get_header().get_tracking().update_time_stamp(t0);
  Field lu(lu_fid);
  lu.allocate_view();
  lu.get_header().get_tracking().update_time_stamp(t0);

  auto time_idx = 0;
  std::vector<Real> iop_tm(nlevs), iop_pm(nlevs), iop_qv(nlevs), iop_qc(nlevs);
  scorpio::read_var(iop_file, "lev", iop_pm.data());
  scorpio::read_var(iop_file, "T", iop_tm.data(), time_idx);
  scorpio::read_var(iop_file, "q", iop_qv.data(), time_idx);
  scorpio::read_var(iop_file, "CLDLIQ", iop_qc.data(), time_idx);
  scorpio::release_file(iop_file);
  scorpio::finalize_subsystem();

  // gravity
  constexpr Real g = PC::gravit;

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0, 0.05);
  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  constexpr int ntests = 1;
  for(int itest = 0; itest < ntests; ++itest) {
    // Randomize
    randomize(qc, engine, pdf);
    randomize(qv, engine, pdf);
    randomize(tm, engine, pdf);
    randomize(pm, engine, pdf);
    randomize(pd, engine, pdf);
    randomize(sd, engine, pdf);
    randomize(su, engine, pdf);
    randomize(ld, engine, pdf);
    randomize(lu, engine, pdf);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    auto diag = diag_factory.create("PBLEntrainmentBudget", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(qc);
    diag->set_required_field(qv);
    diag->set_required_field(tm);
    diag->set_required_field(pm);
    diag->set_required_field(pd);
    diag->set_required_field(sd);
    diag->set_required_field(su);
    diag->set_required_field(ld);
    diag->set_required_field(lu);

    diag->initialize(t0, RunType::Initial);

    qc.get_header().get_tracking().update_time_stamp(t1);
    qv.get_header().get_tracking().update_time_stamp(t1);
    tm.get_header().get_tracking().update_time_stamp(t1);
    pm.get_header().get_tracking().update_time_stamp(t1);
    pd.get_header().get_tracking().update_time_stamp(t1);
    sd.get_header().get_tracking().update_time_stamp(t1);
    su.get_header().get_tracking().update_time_stamp(t1);
    ld.get_header().get_tracking().update_time_stamp(t1);
    lu.get_header().get_tracking().update_time_stamp(t1);

    qc.deep_copy(0.001);
    qv.deep_copy(0.002);
    tm.deep_copy(290.0);
    pm.deep_copy(1e5);
    pd.deep_copy(0.0);
    sd.deep_copy(0.0);
    su.deep_copy(0.0);
    ld.deep_copy(0.0);
    lu.deep_copy(0.0);

    // Run diag
    diag->init_timestep(t1);

    randomize(qc, engine, pdf);
    randomize(qv, engine, pdf);
    randomize(tm, engine, pdf);
    randomize(pm, engine, pdf);
    randomize(pd, engine, pdf);
    randomize(sd, engine, pdf);
    randomize(su, engine, pdf);
    randomize(ld, engine, pdf);
    randomize(lu, engine, pdf);

    // This diag call ensures that random values are
    // correctly set in the diagnostic, with no numerical bugs
    diag->compute_diagnostic();

    auto tm_v = tm.get_view<Real **, Host>();
    auto pm_v = pm.get_view<Real **, Host>();
    auto qc_v = qc.get_view<Real **, Host>();
    auto qv_v = qv.get_view<Real **, Host>();
    auto pd_v = pd.get_view<Real **, Host>();

    for(int ilev = 1; ilev < nlevs; ++ilev) {
      tm_v(0, ilev) = iop_tm[ilev];
      pm_v(0, ilev) = iop_pm[ilev];
      qc_v(0, ilev) = iop_qc[ilev];
      qv_v(0, ilev) = iop_qv[ilev];
      pd_v(0, ilev) = iop_pm[ilev] - iop_pm[ilev - 1];
    }

    tm.sync_to_dev();
    pm.sync_to_dev();
    qc.sync_to_dev();
    qv.sync_to_dev();
    pd.sync_to_dev();

    qc.get_header().get_tracking().update_time_stamp(t2);
    qv.get_header().get_tracking().update_time_stamp(t2);
    tm.get_header().get_tracking().update_time_stamp(t2);
    pm.get_header().get_tracking().update_time_stamp(t2);
    pd.get_header().get_tracking().update_time_stamp(t2);
    sd.get_header().get_tracking().update_time_stamp(t2);
    su.get_header().get_tracking().update_time_stamp(t2);
    ld.get_header().get_tracking().update_time_stamp(t2);
    lu.get_header().get_tracking().update_time_stamp(t2);

    diag->compute_diagnostic();

    // Check result
    qc.sync_to_host();
    qv.sync_to_host();
    tm.sync_to_host();
    pm.sync_to_host();

    diag->get_diagnostic().sync_to_host();

    const auto qc_h = qc.get_view<const Real **, Host>();
    const auto qv_h = qv.get_view<const Real **, Host>();
    const auto tm_h = tm.get_view<const Real **, Host>();
    const auto pm_h = pm.get_view<const Real **, Host>();

    const auto out_hf = diag->get_diagnostic();
    auto out_hv       = out_hf.get_view<Real **, Host>();

    // For the sake of cuteness, let's print out the details, shall we?
    // Set the index map and units map
    PBLEntrainmentBudgetBudgetUtil eadu;
    auto index_map = eadu.index_map;
    auto units_map = eadu.units_map;
    for(const auto &it : index_map) {
      std::cout << it.first << " output is: " << out_hv(0, it.second) << " ("
                << units_map[it.first] << ") \n"
                << std::endl;
    }

    REQUIRE(out_hv(0, 0) == Real(iop_pm[123]));
    std::cout << "Success! Inversion is just below " << iop_pm[123] << " Pa. "
              << "The temperature jump is: " << (iop_tm[123] - iop_tm[124])
              << " K.\n"
              << std::endl;

    REQUIRE(out_hv(0, 4) == Real(iop_qc[123] + iop_qv[123]));
    std::cout << "Success! Calculated qt is " << iop_qc[123] + iop_qv[123]
              << " kg/kg.\n"
              << std::endl;

    Real qt_sum;
    qt_sum = 0.0;
    for(int ilev = 124; ilev < nlevs; ++ilev) {
      qt_sum +=
          (iop_qc[ilev] + iop_qv[ilev]) * (iop_pm[ilev] - iop_pm[ilev - 1]) / g;
    }
    REQUIRE_THAT(out_hv(0, 5), Catch::Matchers::WithinAbs(qt_sum, tol));
    std::cout << "Success! Calculated qt_sum is " << qt_sum << " kg/kg.\n"
              << std::endl;

    auto dt = t2 - t1;
    Real qt_ten;
    qt_ten = 0.0;
    for(int ilev = 124; ilev < nlevs; ++ilev) {
      qt_ten += ((iop_qc[ilev] + iop_qv[ilev] - 0.003) / dt) *
                (iop_pm[ilev] - iop_pm[ilev - 1]) / g;
    }
    REQUIRE_THAT(out_hv(0, 6), Catch::Matchers::WithinAbs(qt_ten, tol));
    std::cout << "Success! Calculated qt_ten is " << qt_ten << " kg/kg/s.\n"
              << std::endl;

  }
}  // TEST_CASE("pbl_entrainment_budget")

}  // namespace scream
