#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {


TEST_CASE("zonal_avg") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using PC = physics::Constants<Real>;

  // A numerical tolerance
  const auto tol = std::numeric_limits<Real>::epsilon() * 100;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  const int nlevs = 3;
  const int dim4  = 4;
  const int ncols = 6;
  const int nlats = 4; // needs to be <= ncols*comm.size()

  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);

  Field area = grid->create_geometry_data("area",grid->get_2d_scalar_layout(),sr);
  area.deep_copy(4*PC::Pi / grid->get_num_global_dofs());
  auto area_view_h = area.get_view<const Real *, Host>();
  area.sync_to_host();

  // Set latitude values
  Field lat = grid->create_geometry_data("lat", grid->get_2d_scalar_layout(), none);
  auto lat_view_h      = lat.get_view<Real *, Host>();
  const Real lat_delta = sp(180.0) / nlats;
  std::vector<Real> zonal_areas(nlats, 0.0);
  for (int i = 0; i < ncols; i++) {
    lat_view_h(i) = sp(-90.0) + (i % nlats + sp(0.5)) * lat_delta;
    zonal_areas[i % nlats] += area_view_h[i];
  }
  comm.all_reduce(zonal_areas.data(), zonal_areas.size(), MPI_SUM);
  lat_view_h(0) = sp(-90.0); // move column to be directly at southern pole
  lat_view_h(nlats-1) = sp(90.0); // move column to be directly at northern pole
  lat.sync_to_dev();

  // Input (randomized) qc
  auto scalar1d_layout = grid->get_2d_scalar_layout();
  auto scalar2d_layout = grid->get_3d_scalar_layout(LEV);
  auto scalar3d_layout = grid->get_3d_vector_layout(LEV,dim4,"dim4");

  FieldIdentifier qc1_fid("qc", scalar1d_layout, kg / kg, grid->name());
  FieldIdentifier qc2_fid("qc", scalar2d_layout, kg / kg, grid->name());
  FieldIdentifier qc3_fid("qc", scalar3d_layout, kg / kg, grid->name());

  Field qc1(qc1_fid);
  Field qc2(qc2_fid);
  Field qc3(qc3_fid);

  qc1.allocate_view();
  qc2.allocate_view();
  qc3.allocate_view();

  // Random number generator seed
  int seed = get_random_test_seed(&comm);

  // Set time for qc and randomize its values
  qc1.get_header().get_tracking().update_time_stamp(t0);
  qc2.get_header().get_tracking().update_time_stamp(t0);
  qc3.get_header().get_tracking().update_time_stamp(t0);
  randomize_uniform(qc1, seed++, 0, 200);
  randomize_uniform(qc2, seed++, 0, 200);
  randomize_uniform(qc3, seed++, 0, 200);

  // Construct the Diagnostics
  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  // Create and set up the diagnostic
  ekat::ParameterList params;
  REQUIRE_THROWS(diag_factory.create("ZonalAvg", comm,
                                     params, grid)); // Bad construction

  params.set("grid_name", grid->name());
  REQUIRE_THROWS(diag_factory.create("ZonalAvg", comm,
                                     params, grid)); // Still no field_name

  params.set<std::string>("field_name", "qc");
  REQUIRE_THROWS(diag_factory.create("ZonalAvg", comm,
                                     params, grid)); // Still no number_of_zonal_bins

  params.set<std::string>("number_of_zonal_bins", std::to_string(nlats));
  // Now we should be good to go...
  auto diag1 = diag_factory.create("ZonalAvg", comm, params, grid);
  auto diag2 = diag_factory.create("ZonalAvg", comm, params, grid);
  auto diag3 = diag_factory.create("ZonalAvg", comm, params, grid);
  // Test the zonal average of qc1
  diag1->set_input_field(qc1);
  diag1->initialize();
  diag1->compute(t0);
  auto diag1_field = diag1->get();

  // Manual calculation
  const std::string bin_dim_name = diag1_field.get_header().get_identifier().get_layout().name(0);
  FieldLayout diag0_layout({CMP}, {nlats}, {bin_dim_name});
  FieldIdentifier diag0_id("qc_zonal_avg_manual", diag0_layout, kg / kg, grid->name());
  Field diag0_field(diag0_id);
  diag0_field.allocate_view();

  // calculate the zonal average
  auto qc1_view_h   = qc1.get_view<const Real *, Host>();
  auto diag0_view_h = diag0_field.get_view<Real *, Host>();
  for (int i = 0; i < ncols; i++) {
    const int nlat = i % nlats;
    diag0_view_h(nlat) += area_view_h(i) / zonal_areas[nlat] * qc1_view_h(i);
  }
  comm.all_reduce(diag0_field.template get_internal_view_data<Real, Host>(),
    diag0_layout.size(), MPI_SUM);
  diag0_field.sync_to_dev();

  // Compare
  REQUIRE(views_are_equal(diag1_field, diag0_field));

  // Try other known cases
  // Set qc1_v to 1.0 to get zonal averages of 1.0/nlats
  const Real zavg1 = sp(1.0);
  qc1.deep_copy(zavg1);
  // Change the input timestamp, to prevent early return and trigger diag recalculation
  t0 += 1;
  qc1.get_header().get_tracking().update_time_stamp(t0);
  diag1->compute(t0);
  auto diag1_view_host = diag1_field.get_view<const Real *, Host>();
  for (int nlat = 0; nlat < nlats; nlat++) {
    REQUIRE_THAT(diag1_view_host(nlat), Catch::Matchers::WithinRel(zavg1, tol));
  }

  // other diags
  // Set qc2_v to 5.0 to get weighted average of 5.0
  const Real zavg2 = sp(5.0);
  qc2.deep_copy(zavg2);
  diag2->set_input_field(qc2);
  diag2->initialize();
  diag2->compute(t0);
  auto diag2_field = diag2->get();

  auto diag2_view_host = diag2_field.get_view<const Real **, Host>();
  for (int i = 0; i < nlevs; ++i) {
    for (int nlat = 0; nlat < nlats; nlat++) {
      REQUIRE_THAT(diag2_view_host(nlat, i), Catch::Matchers::WithinRel(zavg2, tol));
    }
  }

  // Try a random case with qc3
  FieldLayout diag3m_layout({CMP, CMP, LEV}, {nlats, dim4, nlevs},
                            {bin_dim_name, "dim4", e2str(LEV)});
  FieldIdentifier diag3m_id("qc_zonal_avg_manual", diag3m_layout, kg / kg, grid->name());
  Field diag3m_field(diag3m_id);
  diag3m_field.allocate_view();
  auto qc3_view_h    = qc3.get_view<const Real ***, Host>();
  auto diag3m_view_h = diag3m_field.get_view<Real ***, Host>();
  for (int i = 0; i < ncols; i++) {
    const int nlat = i % nlats;
    for (int j = 0; j < dim4; j++) {
      for (int k = 0; k < nlevs; k++) {
        diag3m_view_h(nlat, j, k) += area_view_h(i) / zonal_areas[nlat] * qc3_view_h(i, j, k);
      }
    }
  }
  comm.all_reduce(diag3m_field.template get_internal_view_data<Real, Host>(),
    diag3m_layout.size(), MPI_SUM);
  diag3m_field.sync_to_dev();
  diag3->set_input_field(qc3);
  diag3->initialize();
  diag3->compute(t0);
  auto diag3_field = diag3->get();
  REQUIRE(views_are_equal(diag3_field, diag3m_field));
}

} // namespace scream
