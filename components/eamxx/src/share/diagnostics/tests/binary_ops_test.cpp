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
  FieldIdentifier m_fid("m", scalar2d_layout, kg, grid->name());

  Field qc(qc_fid, true);
  Field qv(qv_fid, true);
  Field m(m_fid, true);

  // Random number generator seed
  int seed = get_random_test_seed();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params_plus, params_times, params_scl_scl;
  REQUIRE_THROWS(diag_factory.create("BinaryOpsDiag", comm,
                                     params_plus));  // No 'arg1', 'arg2', or 'binary_op'

  // Set time for qc and randomize its values
  qc.get_header().get_tracking().update_time_stamp(t0);
  qv.get_header().get_tracking().update_time_stamp(t0);
  m.get_header().get_tracking().update_time_stamp(t0);
  randomize_uniform(qc, seed++, 0, 200);
  randomize_uniform(qv, seed++, 0, 200);
  randomize_uniform(m,  seed++, 0, 200);

  // Create random mask fields for qc/qv
  qc.create_valid_mask();
  qv.create_valid_mask();
  randomize_uniform(qc.get_valid_mask(),seed++);
  randomize_uniform(qv.get_valid_mask(),seed++);

  // Create and set up the diagnostic
  params_plus.set("grid_name", grid->name());
  params_plus.set<std::string>("arg1", "qc");
  params_plus.set<std::string>("arg2", "qv");
  params_plus.set<std::string>("binary_op", "starcraft");
  REQUIRE_THROWS(diag_factory.create("BinaryOpsDiag", comm, params_plus));  // Unknown binary_op
  params_plus.set<std::string>("binary_op", "plus");
  auto plus_diag = diag_factory.create("BinaryOpsDiag", comm, params_plus);

  params_times.set("grid_name", grid->name());
  params_times.set<std::string>("arg1", "m");
  params_times.set<std::string>("arg2", "gravit");
  params_times.set<std::string>("binary_op", "times");
  auto prod_diag = diag_factory.create("BinaryOpsDiag", comm, params_times);

  params_scl_scl.set("grid_name", grid->name());
  params_scl_scl.set<std::string>("arg1", "Avogad");
  params_scl_scl.set<std::string>("arg2", "boltz");
  params_scl_scl.set<std::string>("binary_op", "times");
  auto prod_scl_scl = diag_factory.create("BinaryOpsDiag", comm, params_scl_scl);

  plus_diag->set_grids(gm);
  plus_diag->set_required_field(qc);
  plus_diag->set_required_field(qv);
  plus_diag->initialize(t0, RunType::Initial);

  prod_diag->set_grids(gm);
  prod_diag->set_required_field(m);
  prod_diag->initialize(t0, RunType::Initial);

  prod_scl_scl->set_grids(gm);
  prod_scl_scl->initialize(t0, RunType::Initial);

  // Run diag
  plus_diag->compute_diagnostic();
  prod_diag->compute_diagnostic();
  prod_scl_scl->compute_diagnostic();

  auto plus_diag_f = plus_diag->get_diagnostic(); plus_diag_f.sync_to_host();
  auto prod_diag_f = prod_diag->get_diagnostic(); prod_diag_f.sync_to_host();
  auto rgas_diag_f = prod_scl_scl->get_diagnostic(); rgas_diag_f.sync_to_host();

  // Check that the output fields have the right values
  const auto g = physics::Constants<Real>::gravit.value;

  // field*constant diag
  m.scale(g);
  REQUIRE (not prod_diag_f.has_valid_mask());
  REQUIRE (views_are_equal(prod_diag_f,m));

  // field+field diag
  REQUIRE (plus_diag_f.has_valid_mask());
  auto tgt_mask = qc.get_valid_mask().clone();
  tgt_mask.scale(qv.get_valid_mask());
  REQUIRE (views_are_equal(tgt_mask,plus_diag_f.get_valid_mask()));

  qv.update(qc,1,1,tgt_mask);
  qv.deep_copy(0,tgt_mask,true);
  plus_diag_f.deep_copy(0,tgt_mask,true);
  REQUIRE (views_are_equal(plus_diag_f,qv));

  // constant*constant diag
  REQUIRE (rgas_diag_f.get_view<Real,Host>()()==physics::Constants<Real>::dictionary().at("Rgas").value);

  // redundant, why not
  qc.update(qv, 1, 1);
  views_are_equal(qc, plus_diag_f);
}

}  // namespace scream
