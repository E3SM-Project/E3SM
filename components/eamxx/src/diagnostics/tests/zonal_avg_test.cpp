#include "catch2/catch.hpp"
#include "diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/mesh_free_grids_manager.hpp"
#include "share/util/eamxx_setup_random_test.hpp"
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

TEST_CASE("zonal_avg") {
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
  const int nlats     = 4;

  auto gm   = create_gm(comm, ngcols, nlevs);
  auto grid = gm->get_grid("Physics");

  Field area       = grid->get_geometry_data("area");
  auto area_view_h = area.get_view<const Real *, Host>();

  // Set latitude values
  Field lat = gm->get_grid_nonconst("Physics")->create_geometry_data(
      "lat", grid->get_2d_scalar_layout(), Units::nondimensional());
  auto lat_view_h      = lat.get_view<Real *, Host>();
  const Real lat_delta = sp(180.0) / nlats;
  std::vector<Real> zonal_areas(nlats, 0.0);
  for (int i = 0; i < ngcols; i++) {
    lat_view_h(i) = sp(-90.0) + (i % nlats + sp(0.5)) * lat_delta;
    zonal_areas[i % nlats] += area_view_h[i];
  }
  lat.sync_to_dev();

  // Input (randomized) qc
  FieldLayout scalar1d_layout{{COL}, {ngcols}};
  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
  FieldLayout scalar3d_layout{{COL, CMP, LEV}, {ngcols, dim3, nlevs}};

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
  REQUIRE_THROWS(diag_factory.create("ZonalAvgDiag", comm,
                                     params)); // Bad construction

  params.set("grid_name", grid->name());
  REQUIRE_THROWS(diag_factory.create("ZonalAvgDiag", comm,
                                     params)); // Still no field_name

  params.set<std::string>("field_name", "qc");
  REQUIRE_THROWS(diag_factory.create("ZonalAvgDiag", comm,
                                     params)); // Still no number_of_zonal_bins

  params.set<std::string>("number_of_zonal_bins", std::to_string(nlats));
  // Now we should be good to go...
  auto diag1 = diag_factory.create("ZonalAvgDiag", comm, params);
  auto diag2 = diag_factory.create("ZonalAvgDiag", comm, params);
  auto diag3 = diag_factory.create("ZonalAvgDiag", comm, params);
  diag1->set_grids(gm);
  diag2->set_grids(gm);
  diag3->set_grids(gm);

  // Test the zonal average of qc1
  diag1->set_required_field(qc1);
  diag1->initialize(t0, RunType::Initial);
  diag1->compute_diagnostic();
  auto diag1_field = diag1->get_diagnostic();

  // Manual calculation
  const std::string bin_dim_name = diag1_field.get_header().get_identifier().get_layout().name(0);
  FieldLayout diag0_layout({CMP}, {nlats}, {bin_dim_name});
  FieldIdentifier diag0_id("qc_zonal_avg_manual", diag0_layout, kg / kg, grid->name());
  Field diag0_field(diag0_id);
  diag0_field.allocate_view();

  // calculate the zonal average
  auto qc1_view_h   = qc1.get_view<const Real *, Host>();
  auto diag0_view_h = diag0_field.get_view<Real *, Host>();
  for (int i = 0; i < ngcols; i++) {
    const int nlat = i % nlats;
    diag0_view_h(nlat) += area_view_h(i) / zonal_areas[nlat] * qc1_view_h(i);
  }
  diag0_field.sync_to_dev();

  // Compare
  REQUIRE(views_are_equal(diag1_field, diag0_field));

  // Try other known cases
  // Set qc1_v to 1.0 to get zonal averages of 1.0/nlats
  const Real zavg1 = sp(1.0);
  qc1.deep_copy(zavg1);
  diag1->compute_diagnostic();
  auto diag1_view_host = diag1_field.get_view<Real *, Host>();
  for (int nlat = 0; nlat < nlats; nlat++) {
    REQUIRE_THAT(diag1_view_host(nlat), Catch::Matchers::WithinRel(zavg1, tol));
  }

  // other diags
  // Set qc2_v to 5.0 to get weighted average of 5.0
  const Real zavg2 = sp(5.0);
  qc2.deep_copy(zavg2);
  diag2->set_required_field(qc2);
  diag2->initialize(t0, RunType::Initial);
  diag2->compute_diagnostic();
  auto diag2_field = diag2->get_diagnostic();

  auto diag2_view_host = diag2_field.get_view<Real **, Host>();
  for (int i = 0; i < nlevs; ++i) {
    for (int nlat = 0; nlat < nlats; nlat++) {
      REQUIRE_THAT(diag2_view_host(nlat, i), Catch::Matchers::WithinRel(zavg2, tol));
    }
  }

  // Try a random case with qc3
  FieldLayout diag3m_layout({CMP, CMP, LEV}, {nlats, dim3, nlevs},
                            {bin_dim_name, e2str(CMP), e2str(LEV)});
  FieldIdentifier diag3m_id("qc_zonal_avg_manual", diag3m_layout, kg / kg, grid->name());
  Field diag3m_field(diag3m_id);
  diag3m_field.allocate_view();
  auto qc3_view_h    = qc3.get_view<Real ***, Host>();
  auto diag3m_view_h = diag3m_field.get_view<Real ***, Host>();
  for (int i = 0; i < ngcols; i++) {
    const int nlat = i % nlats;
    for (int j = 0; j < dim3; j++) {
      for (int k = 0; k < nlevs; k++) {
        diag3m_view_h(nlat, j, k) += area_view_h(i) / zonal_areas[nlat] * qc3_view_h(i, j, k);
      }
    }
  }
  diag3m_field.sync_to_dev();
  diag3->set_required_field(qc3);
  diag3->initialize(t0, RunType::Initial);
  diag3->compute_diagnostic();
  auto diag3_field = diag3->get_diagnostic();
  REQUIRE(views_are_equal(diag3_field, diag3m_field));
}

} // namespace scream
