#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_universal_constants.hpp"

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

TEST_CASE("extraaci") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  const auto nondim = Units::nondimensional();

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 5;
  const int ngcols    = 1 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  // Input (randomized) qc, nc
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldIdentifier qc_fid("qc", scalar2d_layout, kg / kg, grid->name());

  Field qc(qc_fid);
  qc.allocate_view();
  qc.get_header().get_tracking().update_time_stamp(t0);

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
      diag_factory.create("AtmTendDiag", comm, params));  // No 'Tend Name'

  //   TODO: The diag currently doesn't throw when given a phony name, need
  //   hardening! params.set<std::string>("Tend Name", "NoWay"); REQUIRE_THROWS
  //   (diag_factory.create("AtmTendDiag",comm,params)); // Bad 'Tend Name'

  // Randomize
  randomize(qc, engine, pdf);
  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("Tend Name", "qc");
  auto diag = diag_factory.create("AtmTendDiag", comm, params);
  diag->set_grids(gm);
  diag->set_required_field(qc);
  diag->initialize(t0, RunType::Initial);

  auto qc_v  = qc.get_view<Real **, Host>();
  qc_v(0, 0) = 5.0;
  qc.sync_to_dev();

  // Run diag
  diag->compute_diagnostic();
  auto diag_f = diag->get_diagnostic();

  Real var_fill_value = constants::DefaultFillValue<Real>().value;

  // Check result: diag should be filled with var_fill_value
  auto some_field = qc.clone();
  some_field.deep_copy(var_fill_value);
  REQUIRE(views_are_equal(diag_f, some_field));

  util::TimeStamp t1({2024, 1, 2}, {0, 0, 0});  // a day later?
  qc.get_header().get_tracking().update_time_stamp(t1);
  qc_v(0, 0) = 29.0;
  qc.sync_to_dev();
  // diag->initialize(t1, RunType::Initial);
  // diag->update(t1, RunType::Initial);
  diag->compute_diagnostic();
  diag_f      = diag->get_diagnostic();
  auto diag_v = diag_f.get_view<Real **, Host>();
  REQUIRE(diag_v(0, 0) == 1.0 / 3600.0);
}

}  // namespace scream
