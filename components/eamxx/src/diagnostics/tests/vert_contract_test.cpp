#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "physics/share/physics_constants.hpp"
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

TEST_CASE("vert_contract") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A numerical tolerance
  auto tol = std::numeric_limits<Real>::epsilon() * 100;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 3;
  constexpr int dim3  = 4;
  const int ngcols    = 6 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) qc
  FieldLayout scalar1d_layout{{LEV}, {nlevs}};
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldLayout scalar3d_layout{{COL, CMP, LEV}, {ngcols, dim3, nlevs}};

  FieldIdentifier fin2_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier fin3_fid("qc", scalar3d_layout, kg / kg, grid->name());
  FieldIdentifier pd_fid("pseudo_density", scalar1d_layout, Pa, grid->name());

  Field fin2(fin2_fid);
  Field fin3(fin3_fid);
  Field pd(pd_fid);

  fin2.allocate_view();
  fin3.allocate_view();
  pd.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(sp(0.0), sp(200.0));
  auto engine = scream::setup_random_test();

  // Construct the diagnostics factory
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("VerContractDiag", comm,
                                     params));  // No 'field_name' parameter

  fin2.get_header().get_tracking().update_time_stamp(t0);
  fin3.get_header().get_tracking().update_time_stamp(t0);
  pd.get_header().get_tracking().update_time_stamp(t0);
  randomize(fin2, engine, pdf);
  randomize(fin3, engine, pdf);
  randomize(pd, engine, pdf);

  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("field_name", "qc");
  params.set<std::string>("contract_method", "avg");
  auto diag_wavg = diag_factory.create("VertContractDiag", comm, params);
  params.set<std::string>("contract_method", "sum");
  auto diag_wsum = diag_factory.create("VertContractDiag", comm, params);
  diag_wavg->set_grids(gm);
  diag_wsum->set_grids(gm);

  using PC         = scream::physics::Constants<Real>;
  constexpr Real g = PC::gravit;
  auto pd_scaled   = pd.clone();
  // scale the area field
  pd_scaled.scale(sp(1.0) / g);
  auto sum = field_sum<Real>(pd_scaled, &comm);
  pd_scaled.scale(sp(1.0) / sum);

  // Fields for manual calculation
  FieldIdentifier diag1_fid("qc_vert_contract_manual",
                            scalar2d_layout.clone().strip_dim(LEV), kg / kg,
                            grid->name());
  FieldIdentifier diag2_fid("qc_vert_contract_manual",
                            scalar3d_layout.clone().strip_dim(LEV), kg / kg,
                            grid->name());

  Field diag1_m(diag1_fid);
  Field diag2_m(diag2_fid);

  diag1_m.allocate_view();
  diag2_m.allocate_view();

  // calculate weighted avg directly
  vert_contraction<Real>(diag1_m, fin2, pd_scaled, &comm);

  // Calculate weighted avg through diags
  diag_wavg->set_required_field(fin2);
  diag_wavg->set_required_field(pd);
  diag_wavg->initialize(t0, RunType::Initial);
  diag_wavg->compute_diagnostic();
  auto diag_wavg_f = diag_wavg->get_diagnostic();

  REQUIRE(views_are_equal(diag_wavg_f, diag1_m));

  // Repeat for another case, but for sum
  auto pd_wsum = pd.clone();
  pd_wsum.scale(sp(1.0) / g);
  vert_contraction<Real>(diag2_m, fin3, pd_wsum, &comm);
  diag_wsum->set_required_field(fin3);
  diag_wsum->set_required_field(pd);
  diag_wsum->initialize(t0, RunType::Initial);
  diag_wsum->compute_diagnostic();
  auto diag_wsum_f = diag_wsum->get_diagnostic();
  REQUIRE(views_are_equal(diag_wsum_f, diag2_m));
}

}  // namespace scream
