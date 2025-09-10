#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include <string>

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ncols, const int nlevs) {
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

TEST_CASE("conditional_sampling") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto fill_value = constants::fill_value<Real>;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 256;
  const int ngcols    = 600 * comm.size();

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("physics");

  // Input (randomized) qc
  FieldLayout scalar1d_layout1{{COL}, {ngcols}};
  FieldLayout scalar1d_layout2{{LEV}, {nlevs}};
  FieldLayout scalar2d_layout1{{COL, LEV}, {ngcols, nlevs}};

  FieldIdentifier qc11_fid("qc", scalar1d_layout1, kg / kg, grid->name());
  FieldIdentifier qc12_fid("qc", scalar1d_layout2, kg / kg, grid->name());
  FieldIdentifier qc21_fid("qc", scalar2d_layout1, kg / kg, grid->name());

  Field qc11(qc11_fid);
  Field qc12(qc12_fid);
  Field qc21(qc21_fid);

  qc11.allocate_view();
  qc12.allocate_view();
  qc21.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(sp(0.0), sp(200.0));
  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  params.set("grid_name", grid->name());

  SECTION("field_v_field") {
    const auto comp_val = 0.001;
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,
                                       params)); // No 'field_name' parameter
    params.set<std::string>("input_field", "qc");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,
                                       params)); // No 'condition_field' parameter
    params.set<std::string>("condition_field", "qc");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,

                                       params)); // No 'condition_operator' parameter
    params.set<std::string>("condition_operator", "gt");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,
                                       params)); // No 'condition_value'
    params.set<std::string>("condition_value", std::to_string(comp_val));

    // Set time for qc and randomize its values
    qc11.get_header().get_tracking().update_time_stamp(t0);
    qc12.get_header().get_tracking().update_time_stamp(t0);
    qc21.get_header().get_tracking().update_time_stamp(t0);
    randomize(qc11, engine, pdf);
    randomize(qc12, engine, pdf);
    randomize(qc21, engine, pdf);

    // Create and set up the diagnostic
    auto diag11 = diag_factory.create("ConditionalSampling", comm, params);
    auto diag12 = diag_factory.create("ConditionalSampling", comm, params);
    auto diag21 = diag_factory.create("ConditionalSampling", comm, params);
    diag11->set_grids(gm);
    diag12->set_grids(gm);
    diag21->set_grids(gm);

    // Set the fields for each diagnostic
    diag11->set_required_field(qc11);
    diag11->initialize(t0, RunType::Initial);
    diag11->compute_diagnostic();
    auto diag11_f = diag11->get_diagnostic();
    diag11_f.sync_to_host();
    auto diag11_v = diag11_f.get_view<const Real *, Host>();

    diag12->set_required_field(qc12);
    diag12->initialize(t0, RunType::Initial);
    diag12->compute_diagnostic();
    auto diag12_f = diag12->get_diagnostic();
    diag12_f.sync_to_host();
    auto diag12_v = diag12_f.get_view<const Real *, Host>();

    diag21->set_required_field(qc21);
    diag21->initialize(t0, RunType::Initial);
    diag21->compute_diagnostic();
    auto diag21_f = diag21->get_diagnostic();
    diag21_f.sync_to_host();
    auto diag21_v = diag21_f.get_view<const Real **, Host>();

    auto qc11_v = qc11.get_view<const Real *, Host>();
    auto qc12_v = qc12.get_view<const Real *, Host>();
    auto qc21_v = qc21.get_view<const Real **, Host>();

    // Check the results
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      // check qc12
      if (qc12_v(ilev) > comp_val) {
        REQUIRE(diag12_v(ilev) == qc12_v(ilev));
      } else {
        REQUIRE(diag12_v(ilev) == fill_value);
      }
    }
    for (int icol = 0; icol < ngcols; ++icol) {
      // Check qc11
      if (qc11_v(icol) > comp_val) {
        REQUIRE(diag11_v(icol) == qc11_v(icol));
      } else {
        REQUIRE(diag11_v(icol) == fill_value);
      }
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        // check qc21
        if (qc21_v(icol, ilev) > comp_val) {
          REQUIRE(diag21_v(icol, ilev) == qc21_v(icol, ilev));
        } else {
          REQUIRE(diag21_v(icol, ilev) == fill_value);
        }
        // check qc21 again, but the negative
        if (qc21_v(icol, ilev) <= comp_val) {
          REQUIRE_FALSE(diag21_v(icol, ilev) == qc21_v(icol, ilev));
        } else {
          REQUIRE_FALSE(diag21_v(icol, ilev) == fill_value);
        }
      }
    }
  }
  SECTION("field_vs_index") {
    const auto comp_lev = static_cast<int>(nlevs / 3);
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,
                                       params)); // No 'field_name' parameter
    params.set<std::string>("input_field", "qc");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,
                                       params)); // No 'condition_field' parameter
    params.set<std::string>("condition_field", "lev");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,

                                       params)); // No 'condition_operator' parameter
    params.set<std::string>("condition_operator", "lt");
    REQUIRE_THROWS(diag_factory.create("ConditionalSampling", comm,
                                       params)); // No 'condition_value'
    params.set<std::string>("condition_value", std::to_string(comp_lev));

    // Set time for qc and randomize its values
    qc11.get_header().get_tracking().update_time_stamp(t0);
    qc12.get_header().get_tracking().update_time_stamp(t0);
    qc21.get_header().get_tracking().update_time_stamp(t0);
    randomize(qc11, engine, pdf);
    randomize(qc12, engine, pdf);
    randomize(qc21, engine, pdf);

    // Create and set up the diagnostic
    auto diag11 = diag_factory.create("ConditionalSampling", comm, params);
    auto diag12 = diag_factory.create("ConditionalSampling", comm, params);
    auto diag21 = diag_factory.create("ConditionalSampling", comm, params);
    diag11->set_grids(gm);
    diag12->set_grids(gm);
    diag21->set_grids(gm);

    // Set the fields for each diagnostic
    diag11->set_required_field(qc11);
    REQUIRE_THROWS(diag11->initialize(t0, RunType::Initial)); // bad dimensions (no lev in qc11)

    diag12->set_required_field(qc12);
    
    diag12->initialize(t0, RunType::Initial);
    diag12->compute_diagnostic();
    auto diag12_f = diag12->get_diagnostic();
    diag12_f.sync_to_host();
    auto diag12_v = diag12_f.get_view<const Real *, Host>();

    diag21->set_required_field(qc21);
    diag21->initialize(t0, RunType::Initial);
    diag21->compute_diagnostic();
    auto diag21_f = diag21->get_diagnostic();
    diag21_f.sync_to_host();
    auto diag21_v = diag21_f.get_view<const Real **, Host>();

    auto qc12_v = qc12.get_view<const Real *, Host>();
    auto qc21_v = qc21.get_view<const Real **, Host>();

    // Check the results
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      // check qc12
      if (ilev < comp_lev) {
        REQUIRE(diag12_v(ilev) == qc12_v(ilev));
      } else {
        REQUIRE(diag12_v(ilev) == fill_value);
      }
    }
    for (int icol = 0; icol < ngcols; ++icol) {
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        // check qc21
        if (ilev < comp_lev) {
          REQUIRE(diag21_v(icol, ilev) == qc21_v(icol, ilev));
        } else {
          REQUIRE(diag21_v(icol, ilev) == fill_value);
        }
        // check qc21 again, but the negative
        if (ilev >= comp_lev) {
          REQUIRE_FALSE(diag21_v(icol, ilev) == qc21_v(icol, ilev));
        } else {
          REQUIRE_FALSE(diag21_v(icol, ilev) == fill_value);
        }
      }
    }
  }
  SECTION("count_conditional") {
    const auto comp_val = 50.0;
    
    // Test count conditional sampling - count grid points where condition is met
    params.set<std::string>("input_field", "count");
    params.set<std::string>("condition_field", "qc");
    params.set<std::string>("condition_operator", "gt");
    params.set<std::string>("condition_value", std::to_string(comp_val));

    // Set time for qc and randomize its values
    qc11.get_header().get_tracking().update_time_stamp(t0);
    qc12.get_header().get_tracking().update_time_stamp(t0);
    qc21.get_header().get_tracking().update_time_stamp(t0);
    randomize(qc11, engine, pdf);
    randomize(qc12, engine, pdf);
    randomize(qc21, engine, pdf);

    // Create and set up the diagnostic for count
    auto count_diag11 = diag_factory.create("ConditionalSampling", comm, params);
    auto count_diag12 = diag_factory.create("ConditionalSampling", comm, params);
    auto count_diag21 = diag_factory.create("ConditionalSampling", comm, params);
    count_diag11->set_grids(gm);
    count_diag12->set_grids(gm);
    count_diag21->set_grids(gm);

    // Set the fields for each diagnostic
    count_diag11->set_required_field(qc11);
    count_diag11->initialize(t0, RunType::Initial);
    count_diag11->compute_diagnostic();
    auto count_diag11_f = count_diag11->get_diagnostic();
    count_diag11_f.sync_to_host();
    auto count_diag11_v = count_diag11_f.get_view<const Real *, Host>();

    count_diag12->set_required_field(qc12);
    count_diag12->initialize(t0, RunType::Initial);
    count_diag12->compute_diagnostic();
    auto count_diag12_f = count_diag12->get_diagnostic();
    count_diag12_f.sync_to_host();
    auto count_diag12_v = count_diag12_f.get_view<const Real *, Host>();

    count_diag21->set_required_field(qc21);
    count_diag21->initialize(t0, RunType::Initial);
    count_diag21->compute_diagnostic();
    auto count_diag21_f = count_diag21->get_diagnostic();
    count_diag21_f.sync_to_host();
    auto count_diag21_v = count_diag21_f.get_view<const Real **, Host>();

    auto qc11_v = qc11.get_view<const Real *, Host>();
    auto qc12_v = qc12.get_view<const Real *, Host>();
    auto qc21_v = qc21.get_view<const Real **, Host>();

    // Check the results - count should be 1.0 where condition is met, 0 otherwise
    for (int ilev = 0; ilev < nlevs; ++ilev) {
      // check count for qc12
      if (qc12_v(ilev) > comp_val) {
        REQUIRE(count_diag12_v(ilev) == 1.0);
      } else {
        REQUIRE(count_diag12_v(ilev) == 0.0);
      }
    }
    
    for (int icol = 0; icol < ngcols; ++icol) {
      // Check count for qc11
      if (qc11_v(icol) > comp_val) {
        REQUIRE(count_diag11_v(icol) == 1.0);
      } else {
        REQUIRE(count_diag11_v(icol) == 0.0);
      }
      
      for (int ilev = 0; ilev < nlevs; ++ilev) {
        // check count for qc21
        if (qc21_v(icol, ilev) > comp_val) {
          REQUIRE(count_diag21_v(icol, ilev) == 1.0);
        } else {
          REQUIRE(count_diag21_v(icol, ilev) == 0.0);
        }
        // check count again, but the negative
        if (qc21_v(icol, ilev) <= comp_val) {
          REQUIRE_FALSE(count_diag21_v(icol, ilev) == 1.0);
        } else {
          REQUIRE_FALSE(count_diag21_v(icol, ilev) == 0.0);
        }
      }
    }
  }
}

} // namespace scream
