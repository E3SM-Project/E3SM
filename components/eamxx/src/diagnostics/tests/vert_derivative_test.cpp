#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols, const int nlevs) {
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

TEST_CASE("vert_derivative") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A numerical tolerance
  // auto tol = std::numeric_limits<Real>::epsilon() * 100;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 150;
  const int ngcols    = 195 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) qc
  FieldLayout scalar1d_layout{{LEV}, {nlevs}};
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};

  FieldIdentifier fin2_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier dp_fid("pseudo_density", scalar2d_layout, Pa, grid->name());

  Field fin2(fin2_fid);
  Field dp(dp_fid);

  fin2.allocate_view();
  dp.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF dp_pdf(10,2000);
  auto engine = scream::setup_random_test();

  // Construct the diagnostics factory
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  // instantiation works because we don't do anything in the constructor
  auto bad_diag = diag_factory.create("VertDerivativeDiag", comm, params);
  SECTION("bad_diag") {
    // this will throw because no field_name was provided
    REQUIRE_THROWS(bad_diag->set_grids(gm));
  }

  fin2.get_header().get_tracking().update_time_stamp(t0);
  dp.get_header().get_tracking().update_time_stamp(t0);
  // instead of random fin2, let's just a predictable one:
  fin2.sync_to_host();
  auto fin2_v_host_assign = fin2.get_view<Real **, Host>();
  for (int icol = 0; icol < ngcols; ++icol) {
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      fin2_v_host_assign(icol, ilev) = 1.0 + icol + ilev;
    }
  }
  fin2.sync_to_dev();
  randomize(dp, engine, dp_pdf);

  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("field_name", "qc");
  SECTION("bad_diag_2") {
    // this will throw because no derivative_method was provided
    auto bad_diag_2 = diag_factory.create("VertDerivativeDiag", comm, params);
    REQUIRE_THROWS(bad_diag_2->set_grids(gm));
  }

  SECTION("bad_diag_3") {
    // this will throw because bad derivative_method was provided
    params.set<std::string>("derivative_method", "xyz");
    auto bad_diag_3 = diag_factory.create("VertDerivativeDiag", comm, params);
    REQUIRE_THROWS(bad_diag_3->set_grids(gm));
  }

  // dp_vert_derivative
  params.set<std::string>("derivative_method", "p");
  auto dp_vert_derivative = diag_factory.create("VertDerivativeDiag", comm, params);

  dp_vert_derivative->set_grids(gm);

  // Fields for manual calculation
  FieldIdentifier diag1_fid("qc_vert_derivative_manual", scalar2d_layout, kg / kg, grid->name());
  Field diag1_m(diag1_fid);
  diag1_m.allocate_view();

  SECTION("dp_vert_derivative") {
    // calculate dp_vert_div manually
    dp.sync_to_host();
    fin2.sync_to_host();
    diag1_m.sync_to_host();

    auto dp_v      = dp.get_view<const Real **, Host>();
    auto fin2_v    = fin2.get_view<const Real **, Host>();
    auto diag1_m_v = diag1_m.get_view<Real **, Host>();

    for (int icol = 0; icol < ngcols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        auto fa1 = (ilev < nlevs - 1) ? (fin2_v(icol, ilev + 1) * dp_v(icol, ilev) +
                    fin2_v(icol, ilev) * dp_v(icol, ilev + 1)) /
                   (dp_v(icol, ilev) + dp_v(icol, ilev + 1)) : fin2_v(icol, nlevs-1);
        auto fa0 = (ilev > 0 ) ? (fin2_v(icol, ilev) * dp_v(icol, ilev - 1) +
                    fin2_v(icol, ilev - 1) * dp_v(icol, ilev)) /
                   (dp_v(icol, ilev - 1) + dp_v(icol, ilev)) : fin2_v(icol, 0);
        diag1_m_v(icol, ilev) = (fa1 - fa0) / dp_v(icol, ilev);
      }
    }
    diag1_m.sync_to_dev();

    // Calculate weighted avg through diagnostics
    dp_vert_derivative->set_required_field(fin2);
    dp_vert_derivative->set_required_field(dp);
    dp_vert_derivative->initialize(t0, RunType::Initial);
    dp_vert_derivative->compute_diagnostic();
    auto dp_vert_derivative_f = dp_vert_derivative->get_diagnostic();

    // Check that the two are the same manually using the tolerance
    dp_vert_derivative_f.sync_to_host();
    auto view_deriv_f_d = dp_vert_derivative_f.get_view<Real **, Host>();
    auto view_deriv_f_m = diag1_m.get_view<Real **, Host>();
    for (int icol = 0; icol < ngcols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        // TODO: the calculations are sometimes diverging because of fp issues in the test on h100 gcc13
        // TODO: so a cop-out, just test if the abs diff is within 1e-5
        REQUIRE(std::abs(view_deriv_f_d(icol, ilev) - view_deriv_f_m(icol, ilev)) < 1e-5);
      }
    }
  }

  // TODO: add SECTION("dz_vert_derivative")
}

} // namespace scream
