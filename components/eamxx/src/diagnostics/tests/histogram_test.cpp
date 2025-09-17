#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

std::shared_ptr<GridsManager> create_gm(const ekat::Comm &comm, const int ngcols, const int nlevs) {
  using vos_t = std::vector<std::string>;
  ekat::ParameterList gm_params;
  gm_params.set("grids_names", vos_t{"Point Grid"});
  auto &pl = gm_params.sublist("Point Grid");
  pl.set<std::string>("type", "point_grid");
  pl.set("aliases", vos_t{"Physics"});
  pl.set<int>("number_of_global_columns", ngcols);
  pl.set<int>("number_of_vertical_levels", nlevs);

  auto gm = create_mesh_free_grids_manager(comm, gm_params);
  gm->build_grids();

  return gm;
}

TEST_CASE("histogram") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A numerical tolerance
  const auto tol = std::numeric_limits<Real>::epsilon() * 100;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  constexpr int nlevs = 3;
  constexpr int dim3  = 4;
  const int ncols     = 6;

  const int ngcols    = ncols * comm.size();
  auto gm             = create_gm(comm, ngcols, nlevs);
  auto grid           = gm->get_grid("Physics");

  // Specify histogram bins
  const std::string bin_configuration = "50_100_150";
  const std::vector<Real> bin_values = {-100.0, 50.0, 100.0, 150.0, 1000.0};

  // Input (randomized) qc
  FieldLayout scalar1d_layout{{COL}, {ncols}};
  FieldLayout scalar2d_layout{{COL, LEV}, {ncols, nlevs}};
  FieldLayout scalar3d_layout{{COL, CMP, LEV}, {ncols, dim3, nlevs}};

  FieldIdentifier qc1_id("qc", scalar1d_layout, kg / kg, grid->name());
  FieldIdentifier qc2_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qc3_fid("qc", scalar3d_layout, kg / kg, grid->name());

  Field qc1(qc1_id);
  Field qc2(qc2_fid);
  Field qc3(qc3_fid);

  qc1.allocate_view();
  qc2.allocate_view();
  qc3.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(sp(0.0), sp(200.0));
  auto engine = scream::setup_random_test();

  // Set time for qc and randomize its values
  qc1.get_header().get_tracking().update_time_stamp(t0);
  qc2.get_header().get_tracking().update_time_stamp(t0);
  qc3.get_header().get_tracking().update_time_stamp(t0);
  randomize(qc1, engine, pdf);
  randomize(qc2, engine, pdf);
  randomize(qc3, engine, pdf);

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  // Create and set up the diagnostic
  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("HistogramDiag", comm,
                                     params)); // Bad construction

  params.set("grid_name", grid->name());
  REQUIRE_THROWS(diag_factory.create("HistogramDiag", comm,
                                     params)); // Still no field_name

  params.set<std::string>("field_name", "qc");
  REQUIRE_THROWS(diag_factory.create("HistogramDiag", comm,
                                     params)); // Still no bin configuration

  params.set<std::string>("bin_configuration", "100_25");
  REQUIRE_THROWS(diag_factory.create("HistogramDiag", comm,
                                     params)); // Non-monotonic bin configuation

  params.set<std::string>("bin_configuration", bin_configuration);

  // Now we should be good to go...
  auto diag1 = diag_factory.create("HistogramDiag", comm, params);
  auto diag2 = diag_factory.create("HistogramDiag", comm, params);
  auto diag3 = diag_factory.create("HistogramDiag", comm, params);
  diag1->set_grids(gm);
  diag2->set_grids(gm);
  diag3->set_grids(gm);

  // Test the zonal average of qc1
  diag1->set_required_field(qc1);
  diag1->initialize(t0, RunType::Initial);
  diag1->compute_diagnostic();
  auto diag1_field = diag1->get_diagnostic();

  // Manual calculation
  const int num_bins = bin_values.size()-1;
  const std::string bin_dim_name = diag1_field.get_header().get_identifier().get_layout().name(0);
  FieldLayout diag0_layout({CMP}, {num_bins}, {bin_dim_name});
  FieldIdentifier diag0_id("qc_histogram_manual", diag0_layout,
    FieldIdentifier::Units::nondimensional(), grid->name());
  Field diag0_field(diag0_id);
  diag0_field.allocate_view();

  // calculate the histogram
  diag0_field.deep_copy(sp(0.0));
  auto qc1_view_h   = qc1.get_view<const Real *, Host>();
  auto diag0_view_h = diag0_field.get_view<Real *, Host>();
  for (int bin_i = 0; bin_i < num_bins; bin_i++) {
    for (int col_i = 0; col_i < ncols; col_i++) {
      if (bin_values[bin_i] <= qc1_view_h(col_i) && qc1_view_h(col_i) < bin_values[bin_i+1])
        diag0_view_h(bin_i) += sp(1.0);
    }
  }
  comm.all_reduce(diag0_field.template get_internal_view_data<Real, Host>(),
    diag0_layout.size(), MPI_SUM);
  diag0_field.sync_to_dev();

  // Compare
  REQUIRE(views_are_equal(diag1_field, diag0_field));

  // Try other known cases
  // Set qc1_v to so histogram is all entries in first bin
  const Real zavg1 = sp(0.5*(bin_values[0]+bin_values[1]));
  qc1.deep_copy(zavg1);
  // Change input timestamp, to prevent early return and trigger diag recalculation
  qc1.get_header().get_tracking().update_time_stamp(t0+1);
  diag1->compute_diagnostic();
  auto diag1_view_host = diag1_field.get_view<const Real *, Host>();
  REQUIRE_THAT(diag1_view_host(0), Catch::Matchers::WithinRel(ngcols, tol));
  for (int bin_i = 1; bin_i < num_bins; bin_i++) {
    REQUIRE(diag1_view_host(bin_i) == sp(0.0));
  }

  // other diags
  // Set qc2_v so histogram is all entries in last bin
  const Real zavg2 = sp(0.5*(bin_values[num_bins-1]+bin_values[num_bins]));
  qc2.deep_copy(zavg2);
  diag2->set_required_field(qc2);
  diag2->initialize(t0, RunType::Initial);
  diag2->compute_diagnostic();
  auto diag2_field = diag2->get_diagnostic();
  auto diag2_view_host = diag2_field.get_view<const Real *, Host>();
  REQUIRE_THAT(diag2_view_host(num_bins-1), Catch::Matchers::WithinRel(ngcols*nlevs, tol));
  for (int bin_i = num_bins-2; bin_i >=0; bin_i--) {
    REQUIRE(diag2_view_host(bin_i) == sp(0.0));
  }

  // Try a random case with qc3
  FieldLayout diag3m_layout({CMP}, {num_bins}, {bin_dim_name});
  FieldIdentifier diag3m_id("qc_zonal_avg_manual", diag3m_layout,
    FieldIdentifier::Units::nondimensional(), grid->name());
  Field diag3m_field(diag3m_id);
  diag3m_field.allocate_view();
  auto qc3_view_h    = qc3.get_view<const Real ***, Host>();
  auto diag3m_view_h = diag3m_field.get_view<Real *, Host>();
  for (int bin_i = 0; bin_i < num_bins; bin_i++) {
    for (int i = 0; i < ncols; i++) {
      for (int j = 0; j < dim3; j++) {
        for (int k = 0; k < nlevs; k++) {
          if (bin_values[bin_i] <= qc3_view_h(i,j,k) && qc3_view_h(i,j,k) < bin_values[bin_i+1])
            diag3m_view_h(bin_i) += sp(1.0);
        }
      }
    }
  }
  comm.all_reduce(diag3m_field.template get_internal_view_data<Real, Host>(),
    diag3m_layout.size(), MPI_SUM);
  diag3m_field.sync_to_dev();
  diag3->set_required_field(qc3);
  diag3->initialize(t0, RunType::Initial);
  diag3->compute_diagnostic();
  auto diag3_field = diag3->get_diagnostic();
  REQUIRE(views_are_equal(diag3_field, diag3m_field));
}

} // namespace scream
