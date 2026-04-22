#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols,
                                        const int nlevs) {
  const int num_global_cols = ncols * comm.size();

  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names", vos_t{"point_grid"});
  auto &pl = gm_params.sublist("point_grid");
  pl.set<std::string>("type", "point_grid");
  pl.set("aliases", vos_t{"physics"});
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

  // Create a grids manager
  constexpr int nlevs = 30;
  constexpr int dim3  = 5;
  const int ngcols    = 95 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("physics");

  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  int seed = get_random_test_seed();

  // Shared field layouts
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldLayout scalar3d_layout{{COL, CMP, LEV}, {ngcols, dim3, nlevs}};

  SECTION("bad_input") {
    ekat::ParameterList params;
    // No field_name -> throws in create_requests
    REQUIRE_THROWS(diag_factory.create("VertContractDiag", comm, params)->set_grids(gm));

    params.set<std::string>("field_name", "qc");
    params.set("grid_name", grid->name());
    // No contract_method -> throws in create_requests
    REQUIRE_THROWS(diag_factory.create("VertContractDiag", comm, params)->set_grids(gm));

    params.set<std::string>("contract_method", "xyz");
    // Invalid contract_method -> throws in create_requests
    REQUIRE_THROWS(diag_factory.create("VertContractDiag", comm, params)->set_grids(gm));
  }

  SECTION("dp_weighted_sum") {
    // Create and randomize a 3D (COL x CMP x LEV) input field
    Field fin(FieldIdentifier("qc", scalar3d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    Field dp(FieldIdentifier("pseudo_density", scalar2d_layout, Pa, grid->name()));
    dp.allocate_view();
    dp.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(dp, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "sum");
    params.set<std::string>("weighting_method", "dp");
    auto diag = diag_factory.create("VertContractDiag", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(fin);
    diag->set_required_field(dp);
    diag->initialize(t0, RunType::Initial);
    diag->compute_diagnostic();
    auto diag_f = diag->get_diagnostic();

    // Manual reference: ref(col,cmp) = sum_lev( (dp/g)(col,lev) * f(col,cmp,lev) )
    auto dp_scaled = dp.clone("dp_scaled");
    dp_scaled.scale(1 / physics::Constants<Real>::gravit.value);
    auto ref = diag_f.clone("ref");
    vert_contraction(ref, fin, dp_scaled);

    REQUIRE(views_are_equal(diag_f, ref));
  }

  SECTION("dp_weighted_avg") {
    // Create and randomize a 2D (COL x LEV) input field
    Field fin(FieldIdentifier("qc", scalar2d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    Field dp(FieldIdentifier("pseudo_density", scalar2d_layout, Pa, grid->name()));
    dp.allocate_view();
    dp.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(dp, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "avg");
    params.set<std::string>("weighting_method", "dp");
    auto diag = diag_factory.create("VertContractDiag", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(fin);
    diag->set_required_field(dp);
    diag->initialize(t0, RunType::Initial);
    diag->compute_diagnostic();
    auto diag_f = diag->get_diagnostic();

    // Manual reference:
    //   num(col) = sum_lev( (dp/g)(col,lev) * f(col,lev) )
    //   den(col) = sum_lev( (dp/g)(col,lev) )
    //   ref(col) = num(col) / den(col)
    auto dp_scaled = dp.clone("dp_scaled");
    dp_scaled.scale(1 / physics::Constants<Real>::gravit.value);
    auto dp_ones = dp_scaled.clone("ones");
    dp_ones.deep_copy(1);
    auto ref = diag_f.clone("ref");
    auto den = diag_f.clone("den");
    vert_contraction(ref, fin, dp_scaled);
    vert_contraction(den, dp_scaled, dp_ones);
    ref.scale_inv(den);

    REQUIRE(views_are_equal(diag_f, ref));
  }

  SECTION("unweighted_sum") {
    // Create and randomize a 2D (COL x LEV) input field
    Field fin(FieldIdentifier("qc", scalar2d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "sum");
    params.set<std::string>("weighting_method", "none");
    auto diag = diag_factory.create("VertContractDiag", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(fin);
    diag->initialize(t0, RunType::Initial);
    diag->compute_diagnostic();
    auto diag_f = diag->get_diagnostic();

    // Manual reference: ref(col) = sum_lev( f(col,lev) )
    // Use a 1D all-ones weight (same as the diagnostic)
    FieldLayout ones_lev_layout{{LEV}, {nlevs}};
    Field ones_lev(FieldIdentifier("ones_lev", ones_lev_layout,
                                   Units::nondimensional(), grid->name()), true);
    ones_lev.deep_copy(1);
    auto ref = diag_f.clone("ref");
    vert_contraction(ref, fin, ones_lev);

    REQUIRE(views_are_equal(diag_f, ref));
  }

  SECTION("unweighted_avg") {
    // Create and randomize a 3D (COL x CMP x LEV) input field
    Field fin(FieldIdentifier("qc", scalar3d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "avg");
    params.set<std::string>("weighting_method", "none");
    auto diag = diag_factory.create("VertContractDiag", comm, params);
    diag->set_grids(gm);
    diag->set_required_field(fin);
    diag->initialize(t0, RunType::Initial);
    diag->compute_diagnostic();
    auto diag_f = diag->get_diagnostic();

    // Manual reference:
    //   num(col,cmp) = sum_lev( f(col,cmp,lev) )
    //   den(col,cmp) = sum_lev( 1 ) = nlevs
    //   ref(col,cmp) = num(col,cmp) / den(col,cmp)
    // Use a 1D all-ones weight (same as the diagnostic)
    FieldLayout ones_lev_layout{{LEV}, {nlevs}};
    Field ones_lev(FieldIdentifier("ones_lev", ones_lev_layout,
                                   Units::nondimensional(), grid->name()), true);
    ones_lev.deep_copy(1);
    auto fin_ones = fin.clone("fin_ones");
    fin_ones.deep_copy(1);
    auto ref = diag_f.clone("ref");
    auto den = diag_f.clone("den");
    vert_contraction(ref, fin, ones_lev);
    vert_contraction(den, fin_ones, ones_lev);
    ref.scale_inv(den);

    REQUIRE(views_are_equal(diag_f, ref));
  }
}

}  // namespace scream
