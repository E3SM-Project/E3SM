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

TEST_CASE("zonal_avg") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using TeamPolicy = Kokkos::TeamPolicy<Field::device_t::execution_space>;
  using TeamMember = typename TeamPolicy::member_type;
  using KT         = ekat::KokkosTypes<DefaultDevice>;
  using ESU        = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

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

  // Set latitude values
  Field lat = gm->get_grid_nonconst("Physics")->create_geometry_data("lat",
    grid->get_2d_scalar_layout(), ekat::units::Units::nondimensional());
  auto lat_view = lat.get_view<Real *>();
  for (int i=0; i < ngcols; i++)
  {
    lat_view(i) = (i % 2 == 0) ? -45.0 : 45.0;
  }

  // Input (randomized) qc
  FieldLayout scalar1d_layout{{COL}, {ngcols}};
//  FieldLayout scalar2d_layout{{COL, LEV}, {ngcols, nlevs}};
//  FieldLayout scalar3d_layout{{COL, CMP, LEV}, {ngcols, dim3, nlevs}};

  FieldIdentifier qc1_fid("qc", scalar1d_layout, kg / kg, grid->name());
//  FieldIdentifier qc2_fid("qc", scalar2d_layout, kg / kg, grid->name());
//  FieldIdentifier qc3_fid("qc", scalar3d_layout, kg / kg, grid->name());

  Field qc1(qc1_fid);
//  Field qc2(qc2_fid);
//  Field qc3(qc3_fid);

  qc1.allocate_view();
//  qc2.allocate_view();
//  qc3.allocate_view();

  // Construct random number generator stuff
  using RPDF = std::uniform_real_distribution<Real>;
  RPDF pdf(sp(0.0), sp(200.0));
  auto engine = scream::setup_random_test();

  // Construct the Diagnostics
  std::map<std::string, std::shared_ptr<AtmosphereDiagnostic>> diags;
  auto &diag_factory = AtmosphereDiagnosticFactory::instance();
  register_diagnostics();

  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("ZonalAvgDiag", comm,
                                     params));  // No 'field_name' parameter

  // Set time for qc and randomize its values
  qc1.get_header().get_tracking().update_time_stamp(t0);
//  qc2.get_header().get_tracking().update_time_stamp(t0);
//  qc3.get_header().get_tracking().update_time_stamp(t0);
  randomize(qc1, engine, pdf);
//  randomize(qc2, engine, pdf);
//  randomize(qc3, engine, pdf);

  // Create and set up the diagnostic
  params.set("grid_name", grid->name());
  params.set<std::string>("field_name", "qc");
  auto diag1 = diag_factory.create("ZonalAvgDiag", comm, params);
//  auto diag2 = diag_factory.create("ZonalAvgDiag", comm, params);
//  auto diag3 = diag_factory.create("ZonalAvgDiag", comm, params);
  diag1->set_grids(gm);
//  diag2->set_grids(gm);
//  diag3->set_grids(gm);

  // Clone the area field
  Field scaled_area = grid->get_geometry_data("area").clone();

  // Test the zonal average of qc1
  diag1->set_required_field(qc1);
  diag1->initialize(t0, RunType::Initial);
  diag1->compute_diagnostic();
  auto diag1_field = diag1->get_diagnostic();

  // Manual calculation
  FieldLayout diag0_layout({COL}, {2}, {"lat"});
  FieldIdentifier diag0_id("qc_zonal_avg_manual", diag0_layout, kg / kg,
    grid->name());
  Field diag0_field(diag0_id);
  diag0_field.allocate_view();

  // calculate total area
  Real atot = field_sum<Real>(scaled_area, &comm);
  // scale the area field
  scaled_area.scale(1.0 / atot);

  // calculate the zonal average
  auto qc1_view = qc1.get_view<const Real *>();
  auto diag0_view = diag0_field.get_view<Real *>();
  auto scaled_area_view = scaled_area.get_view<const Real *>();
  for (int i=0; i < ngcols; i++)
  {
    diag0_view(i%2) += scaled_area_view(i) * qc1_view(i);
  }

  // Compare
  REQUIRE(views_are_equal(diag1_field, diag0_field));

  // Try other known cases
  // Set qc1_v to 1.0 to get zonal averages of 0.5
  Real wavg = 1.0;
  qc1.deep_copy(wavg);
  diag1->compute_diagnostic();
  auto diag1_v2_host = diag1_field.get_view<Real *, Host>();
  REQUIRE_THAT(diag1_v2_host(0), Catch::Matchers::WithinRel(0.5*wavg, tol));
  REQUIRE_THAT(diag1_v2_host(1), Catch::Matchers::WithinRel(0.5*wavg, tol));
/*
  // other diags
  // Set qc2_v to 5.0 to get weighted average of 5.0
  wavg = sp(5.0);
  qc2.deep_copy(wavg);
  diag2->set_required_field(qc2);
  diag2->initialize(t0, RunType::Initial);
  diag2->compute_diagnostic();
  auto diag2_f = diag2->get_diagnostic();

  auto diag2_v_host = diag2_f.get_view<Real *, Host>();

  for(int i = 0; i < nlevs; ++i) {
    REQUIRE_THAT(diag2_v_host(i), Catch::Matchers::WithinRel(wavg, tol));
  }

  // Try a random case with qc3
  auto qc3_v = qc3.get_view<Real ***>();
  FieldIdentifier diag3_manual_fid("qc_zonal_avg_manual",
                                   scalar3d_layout.clone().strip_dim(COL),
                                   kg / kg, grid->name());
  Field diag3_manual(diag3_manual_fid);
  diag3_manual.allocate_view();
  horiz_contraction<Real>(diag3_manual, qc3, area, &comm);
  diag3->set_required_field(qc3);
  diag3->initialize(t0, RunType::Initial);
  diag3->compute_diagnostic();
  auto diag3_f = diag3->get_diagnostic();
  REQUIRE(views_are_equal(diag3_f, diag3_manual));
*/
}

}  // namespace scream
