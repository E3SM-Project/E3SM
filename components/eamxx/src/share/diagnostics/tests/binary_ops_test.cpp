#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/data_managers/library_grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols,
                                        const int nlevs) {
  const int num_global_cols = ncols * comm.size();

  auto grid = create_point_grid("physics",num_global_cols,nlevs,comm);
  auto gm = std::make_shared<LibraryGridsManager>();
  gm->add_grids(grid);
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
  const auto &plus_v = plus_diag_f.get_view<Real**, Host>();
  const auto &prod_v = prod_diag_f.get_view<Real**, Host>();
  const auto &rgas_v = rgas_diag_f.get_view<Real, Host>();
  const auto &qc_v   = qc.get_view<Real**, Host>();
  const auto &qv_v   = qv.get_view<Real**, Host>();
  const auto &m_v    = m.get_view<Real**, Host>();
  const auto g = physics::Constants<Real>::gravit.value;
  for (int icol = 0; icol < ngcols; ++icol) {
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      // Check plus
      REQUIRE(plus_v(icol, ilev) == qc_v(icol, ilev) + qv_v(icol, ilev));
      // Check product
      REQUIRE(prod_v(icol, ilev) == (m_v(icol,ilev)*g));
    }
  }
  // The diag_scl_scl shoould compute the prod of avogadro's number and boltzman's constant, which is Rgas
  REQUIRE (rgas_v()==physics::Constants<Real>::dictionary().at("Rgas").value);

  // redundant, why not
  qc.update(qv, 1, 1);
  views_are_equal(qc, plus_diag_f);
}

TEST_CASE ("inputs_have_mask") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});
  const int ncols = 5;

  int seed = get_random_test_seed();

  auto gm = create_gm(comm,ncols,1);
  auto grid = gm->get_grid("physics");

  auto layout = grid->get_2d_scalar_layout();
  FieldIdentifier fid1("f1",layout,m/s,grid->name());
  FieldIdentifier fid2("f2",layout,m/s,grid->name());
  FieldIdentifier mfid("m",layout,m/s,grid->name(),DataType::IntType);

  auto create_diag = [&](const Field& f1, const Field& f2) {
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("arg1", "f1");
    params.set<std::string>("arg2", "f2");
    params.set<std::string>("binary_op", "times");

    auto diag = std::make_shared<BinaryOpsDiag>(comm, params);
    diag->set_grids(gm);

    diag->set_required_field(f1);
    diag->set_required_field(f2);

    diag->initialize(t0,RunType::Initial);
    return diag;
  };

  constexpr auto fv = constants::fill_value<Real>;

  typename Field::view_host_t<const int*> empty;
  for (auto mask1 : {false, true}) {
    Field f1 (fid1,true);
    f1.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(f1,seed++,0,1);

    Field m1;
    if (mask1) {
      m1 = Field(mfid,true);
      randomize_uniform(m1,seed++,0,1);
      f1.get_header().set_extra_data("valid_mask",m1);
    }

    for (auto mask2 : {false, true}) {
      Field f2 (fid2,true);
      f2.get_header().get_tracking().update_time_stamp(t0);
      randomize_uniform(f2,seed++,0,1);

      Field m2;
      if (mask2) {
        m2 = Field(mfid,true);
        randomize_uniform(m2,seed++,0,1);
        f2.get_header().set_extra_data("valid_mask",m2);
      }

      auto diag = create_diag(f1,f2);
      diag->compute_diagnostic();

      auto d = diag->get_diagnostic();
      d.sync_to_host();
      auto dm = mask1 or mask2 ? d.get_header().get_extra_data<Field>("valid_mask") : Field{};

      auto dh  = d.get_view<const Real*,Host>();
      auto f1h = f1.get_view<const Real*,Host>();
      auto f2h = f2.get_view<const Real*,Host>();

      auto dmh = mask1 or mask2 ? dm.get_view<const int*,Host>() : empty;
      auto m1h = mask1 ? m1.get_view<const int*,Host>() : empty;
      auto m2h = mask2 ? m2.get_view<const int*,Host>() : empty;

      if (mask1 and mask2) {
        for (int i=0; i<ncols; ++i) {
          if (m1h[i] and m2h[i])
            REQUIRE (dh[i]==f1h[i]*f2h[i]);
          else
            REQUIRE (dh[i]==fv);

          // Also check that diag mask is the prod of input masks
          REQUIRE (dmh[i]==m1h[i]*m2h[i]);
        }
      } else if (mask1) {
        for (int i=0; i<ncols; ++i) {
          if (m1h[i])
            REQUIRE (dh[i]==f1h[i]*f2h[i]);
          else
            REQUIRE (dh[i]==fv);
        }
        // Diag mask should alias f1's
        REQUIRE (dm.is_aliasing(m1));
      } else if (mask2) {
        for (int i=0; i<ncols; ++i) {
          if (m2h[i])
            REQUIRE (dh[i]==f1h[i]*f2h[i]);
          else
            REQUIRE (dh[i]==fv);
        }
        // Diag mask should alias f2's
        REQUIRE (dm.is_aliasing(m2));
      } else {
        for (int i=0; i<ncols; ++i) {
          REQUIRE (dh[i]==f1h[i]*f2h[i]);
        }
      }
    }
  }
}

}  // namespace scream
