#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
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

TEST_CASE("vert_contract") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 30;
  constexpr int dim3  = 5;
  const int ngcols    = 95 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) qc
  FieldLayout scalar1d_layout{{LEV}, {nlevs}};
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldLayout scalar3d_layout{{COL, CMP, LEV}, {ngcols, dim3, nlevs}};

  FieldIdentifier fin2_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier fin3_fid("qc", scalar3d_layout, kg / kg, grid->name());
  FieldIdentifier dp_fid("pseudo_density", scalar2d_layout, Pa, grid->name());
  FieldIdentifier dz_fid("dz", scalar2d_layout, m, grid->name());

  Field fin2(fin2_fid);
  Field fin3(fin3_fid);
  Field dp(dp_fid);
  // Field dz(dz_fid);

  fin2.allocate_view();
  fin3.allocate_view();
  dp.allocate_view();
  // dz.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(sp(0.0), sp(200.0));
  auto engine = scream::setup_random_test();

  // Construct the diagnostics factory
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  // instantiation works because we don't do anything in the constructor
  auto bad_diag = diag_factory.create("VertContractDiag", comm, params);
  SECTION("bad_diag") {
    // this will throw because no field_name was provided
    REQUIRE_THROWS(bad_diag->set_grids(gm));
  }

  fin2.get_header().get_tracking().update_time_stamp(t0);
  fin3.get_header().get_tracking().update_time_stamp(t0);
  dp.get_header().get_tracking().update_time_stamp(t0);
  // dz.get_header().get_tracking().update_time_stamp(t0);
  randomize(fin2, engine, pdf);
  randomize(fin3, engine, pdf);
  randomize(dp, engine, pdf);
  // randomize(dz, engine, pdf);

  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("field_name", "qc");
  SECTION("bad_diag_2") {
    // this will throw because no contract_method was provided
    auto bad_diag_2 = diag_factory.create("VertContractDiag", comm, params);
    REQUIRE_THROWS(bad_diag_2->set_grids(gm));
  }

  SECTION("bad_diag_3") {
    // this will throw because bad contract_method was provided
    params.set<std::string>("contract_method", "xyz");
    auto bad_diag_3 = diag_factory.create("VertContractDiag", comm, params);
    REQUIRE_THROWS(bad_diag_3->set_grids(gm));
  }

  // dp_weighted_avg
  params.set<std::string>("contract_method", "avg");
  params.set<std::string>("weighting_method", "dp");
  auto dp_weighted_avg = diag_factory.create("VertContractDiag", comm, params);

  // dp_weighted_sum
  params.set<std::string>("contract_method", "sum");
  params.set<std::string>("weighting_method", "dp");
  auto dp_weighted_sum = diag_factory.create("VertContractDiag", comm, params);

  // TODO: for some reason the dz field keeps getting set to 0
  // TODO: as a workaround, just calculate dz here (sigh...)
  // TODO: commenting out all the dz stuff for now
  // // dz_weighted_avg
  // params.set<std::string>("contract_method", "avg");
  // params.set<std::string>("weighting_method", "dz");
  // auto dz_weighted_avg = diag_factory.create("VertContractDiag", comm, params);

  // // dz_weighted_sum
  // params.set<std::string>("contract_method", "sum");
  // params.set<std::string>("weighting_method", "dz");
  // auto dz_weighted_sum = diag_factory.create("VertContractDiag", comm, params);
  
  // unweighted_sum
  params.set<std::string>("contract_method", "sum");
  params.set<std::string>("weighting_method", "none");
  auto unweighted_sum = diag_factory.create("VertContractDiag", comm, params);
  
  // unweighted_avg
  params.set<std::string>("contract_method", "avg");
  params.set<std::string>("weighting_method", "none");
  auto unweighted_avg = diag_factory.create("VertContractDiag", comm, params);

  dp_weighted_avg->set_grids(gm);
  dp_weighted_sum->set_grids(gm);
  // dz_weighted_sum->set_grids(gm);
  // dz_weighted_avg->set_grids(gm);
  unweighted_sum->set_grids(gm);
  unweighted_avg->set_grids(gm);

  // Fields for manual calculation
  FieldIdentifier diag1_fid("qc_vert_contract_manual", scalar2d_layout.clone().strip_dim(LEV), kg / kg, grid->name());
  FieldIdentifier diag2_fid("qc_vert_contract_manual", scalar3d_layout.clone().strip_dim(LEV), kg / kg, grid->name());
  Field diag1_m(diag1_fid);
  Field diag2_m(diag2_fid);
  diag1_m.allocate_view();
  diag2_m.allocate_view();

  // Fields for scaling
  FieldIdentifier dps_fid ("dps", scalar2d_layout.clone().strip_dim(LEV), Pa, grid->name());
  // FieldIdentifier dzs_fid ("dzs", scalar2d_layout.clone().strip_dim(LEV), m, grid->name());
  Field dps(dps_fid);
  // Field dzs(dzs_fid);
  dps.allocate_view();
  // dzs.allocate_view();

  auto dp_ones = dp.clone("dp_ones");
  dp_ones.deep_copy(1);
  // auto dz_ones = dz.clone("dz_ones");
  // dz_ones.deep_copy(1);

  auto dp_scaled   = dp.clone("dp_scaled");
  // auto dz_scaled   = dz.clone("dz_scaled");

  dp_scaled.scale(sp(1.0) / scream::physics::Constants<Real>::gravit);

  vert_contraction<Real>(dps, dp_scaled, dp_ones);
  // vert_contraction<Real>(dzs, dz_scaled, dz_ones);

  SECTION("dp_weighted_avg") {
    // scale dp_scaled by 1/dps (because we are averaging)
    dps.sync_to_host();
    auto dps_v = dps.get_view<const Real*, Host>();
    dp_scaled.sync_to_host();
    auto dp_scaled_v = dp_scaled.get_view<Real**, Host>();
    for (std::size_t i = 0; i < dp_scaled_v.extent(0); ++i) {
      for (std::size_t j = 0; j < dp_scaled_v.extent(1); ++j) {
        if(dps_v(i) == 0) {
          dp_scaled_v(i, j) = 0; // Handle division by zero by setting to 0
        } else {
          dp_scaled_v(i, j) /= dps_v(i);
        }
      }
    }
    dp_scaled.sync_to_dev();

    // calculate weighted avg directly
    vert_contraction<Real>(diag1_m, fin2, dp_scaled);

    // Calculate weighted avg through diagnostics
    dp_weighted_avg->set_required_field(fin2);
    dp_weighted_avg->set_required_field(dp);
    dp_weighted_avg->initialize(t0, RunType::Initial);
    dp_weighted_avg->compute_diagnostic();
    auto dp_weighted_avg_f = dp_weighted_avg->get_diagnostic();

    REQUIRE(views_are_equal(dp_weighted_avg_f, diag1_m));
  }

  SECTION("dp_weighted_sum") {
    // calculate weighted sum directly
    vert_contraction<Real>(diag2_m, fin3, dp_scaled);
    // Calculate weighted sum through diagnostics
    dp_weighted_sum->set_required_field(fin3);
    dp_weighted_sum->set_required_field(dp);
    dp_weighted_sum->initialize(t0, RunType::Initial);
    dp_weighted_sum->compute_diagnostic();
    auto dp_weighted_sum_f = dp_weighted_sum->get_diagnostic();

    REQUIRE(views_are_equal(dp_weighted_sum_f, diag2_m));
  }

  // SECTION("dz_weighted_avg") {
  //   // scale dz_scaled by 1/dzs (because we are averaging)
  //   dzs.sync_to_host();
  //   auto dzs_v = dzs.get_view<const Real*, Host>();
  //   dz_scaled.sync_to_host();
  //   auto dz_scaled_v = dz_scaled.get_view<Real**, Host>();
  //   for (std::size_t i = 0; i < dz_scaled_v.extent(0); ++i) {
  //     for (std::size_t j = 0; j < dz_scaled_v.extent(1); ++j) {
  //       if(dzs_v(i) == 0) {
  //         dz_scaled_v(i, j) = 0; // Handle division by zero by setting to 0
  //       } else {
  //         dz_scaled_v(i, j) /= dzs_v(i);
  //       }
  //     }
  //   }
  //   dz_scaled.sync_to_dev();

  //   // calculate weighted avg directly
  //   vert_contraction<Real>(diag1_m, fin2, dz_scaled);

  //   // Calculate weighted avg through diagnostics
  //   dz_weighted_avg->set_required_field(fin2);
  //   dz_weighted_avg->set_required_field(dz);
  //   dz_weighted_avg->initialize(t0, RunType::Initial);
  //   dz_weighted_avg->compute_diagnostic();
  //   auto dz_weighted_avg_f = dz_weighted_avg->get_diagnostic();

  //   REQUIRE(views_are_equal(dz_weighted_avg_f, diag1_m));
  // }

  // SECTION("dz_weighted_sum") {
  //   // calculate weighted sum directly
  //   vert_contraction<Real>(diag2_m, fin3, dz_scaled);
  //   // Calculate weighted sum through diagnostics
  //   dz_weighted_sum->set_required_field(fin3);
  //   dz_weighted_sum->set_required_field(dz);
  //   dz_weighted_sum->initialize(t0, RunType::Initial);
  //   dz_weighted_sum->compute_diagnostic();
  //   auto dz_weighted_sum_f = dz_weighted_sum->get_diagnostic();

  //   REQUIRE(views_are_equal(dz_weighted_sum_f, diag2_m));
  // }

  SECTION("unweighted_sum") {
    // calculate unweighted sum directly
    vert_contraction<Real>(diag1_m, fin2, dp_ones);

    // Calculate unweighted sum through diagnostics
    unweighted_sum->set_required_field(fin2);
    unweighted_sum->initialize(t0, RunType::Initial);
    unweighted_sum->compute_diagnostic();
    auto unweighted_sum_f = unweighted_sum->get_diagnostic();

    REQUIRE(views_are_equal(unweighted_sum_f, diag1_m));
  }

  SECTION("unweighted_avg") {
    // since we are averaging, we need to scale by the sum
    auto dp_ones_scaled = dp_ones.clone("dz_ones_scaled");
    dp_ones_scaled.sync_to_host();
    auto dp_ones_scaled_v = dp_ones_scaled.get_view<Real**, Host>();
    for (std::size_t i = 0; i < dp_ones_scaled_v.extent(0); ++i) {
      for (std::size_t j = 0; j < dp_ones_scaled_v.extent(1); ++j) {
        const int nlevs = dp_ones_scaled_v.extent(1);  
        dp_ones_scaled_v(i, j) /= nlevs;
      }
    }
    dp_ones_scaled.sync_to_dev();
    // calculate unweighted avg directly
    vert_contraction<Real>(diag2_m, fin3, dp_ones_scaled);
    // Calculate unweighted avg through diagnostics
    unweighted_avg->set_required_field(fin3);
    unweighted_avg->initialize(t0, RunType::Initial);
    unweighted_avg->compute_diagnostic();
    auto unweighted_avg_f = unweighted_avg->get_diagnostic();
    
    REQUIRE(views_are_equal(unweighted_avg_f, diag2_m));
  }
}

}  // namespace scream
