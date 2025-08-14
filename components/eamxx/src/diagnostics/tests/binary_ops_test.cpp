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

TEST_CASE("binary_ops") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 201;
  const int ngcols    = 260 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("physics");

  // Input (randomized) qc, qv
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qv_fid("qv", scalar2d_layout, kg / kg, grid->name());

  Field qc(qc_fid);
  Field qv(qv_fid);
  qc.allocate_view();
  qv.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(0.0, 200.0);

  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("BinaryOpsDiag", comm,
                                     params));  // No 'field_1', 'field_2', or 'binary_op'

  // Set time for qc and randomize its values
  qc.get_header().get_tracking().update_time_stamp(t0);
  qv.get_header().get_tracking().update_time_stamp(t0);
  randomize(qc, engine, pdf); qc.sync_to_dev();
  randomize(qv, engine, pdf); qv.sync_to_dev();

  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("field_1", "qc");
  params.set<std::string>("field_2", "qv");
  params.set<std::string>("binary_op", "plus");
  auto plus_diag = diag_factory.create("BinaryOpsDiag", comm, params);
  params.set<std::string>("binary_op", "times");
  auto prod_diag = diag_factory.create("BinaryOpsDiag", comm, params);
  plus_diag->set_grids(gm);
  prod_diag->set_grids(gm);
  plus_diag->set_required_field(qc);
  prod_diag->set_required_field(qc);
  plus_diag->set_required_field(qv);
  prod_diag->set_required_field(qv);
  plus_diag->initialize(t0, RunType::Initial);
  prod_diag->initialize(t0, RunType::Initial);

  // Run diag
  plus_diag->compute_diagnostic();
  auto plus_diag_f = plus_diag->get_diagnostic(); plus_diag_f.sync_to_host();
  prod_diag->compute_diagnostic();
  auto prod_diag_f = prod_diag->get_diagnostic(); prod_diag_f.sync_to_host();

  // Check that the output fields have the right values
  const auto &plus_v = plus_diag_f.get_view<Real**, Host>();
  const auto &prod_v = prod_diag_f.get_view<Real**, Host>();
  const auto &qc_v   = qc.get_view<Real**, Host>();
  const auto &qv_v   = qv.get_view<Real**, Host>();
  for (int icol = 0; icol < ngcols; ++icol) {
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      // Check plus
      REQUIRE(plus_v(icol, ilev) == qc_v(icol, ilev) + qv_v(icol, ilev));
      // Check product
      REQUIRE(prod_v(icol, ilev) == qc_v(icol, ilev) * qv_v(icol, ilev));
    }
  }

  // redundant, why not
  qc.update(qv, 1, 1);
  views_are_equal(qc, plus_diag_f);
}

}  // namespace scream
