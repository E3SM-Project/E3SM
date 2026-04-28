#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {


TEST_CASE("horiz_avg") {
  using namespace ShortFieldTagsNames;
  using PC = physics::Constants<Real>;

  // A numerical tolerance
  auto tol = std::numeric_limits<Real>::epsilon() * 100;

  // Random number generator seed
  int seed = get_random_test_seed();

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager - single column for these tests
  const int nlevs = 3;
  const int ncols = 6;

  auto grid = create_point_grid("physics",ncols,nlevs,comm);
  auto area = grid->create_geometry_data("area",grid->get_2d_scalar_layout());
  area.deep_copy(4.0*PC::Pi/grid->get_num_global_dofs());

  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  // No 'field_name' parameter
  REQUIRE_THROWS(diag_factory.create("HorizAvg", comm, ekat::ParameterList(), grid));

  // Get a hold to area, for easy cloning later
  auto area_sum = field_sum(area,&comm).as<Real>();

  SECTION ("scalar_x") {
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "s_x");

    FieldLayout layout = grid->get_2d_scalar_layout();
    FieldIdentifier scl_x_fid ("s_x" , layout, ekat::units::none, grid->name());
    Field scl_x (scl_x_fid , true);
    scl_x.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(scl_x,seed++);

    auto diag = diag_factory.create("HorizAvg", comm, params, grid);
    diag->set_input_field(scl_x);
    diag->initialize();
    diag->compute_diagnostic(t0);
    auto diag_f = diag->get_diagnostic();

    auto manual = diag_f.clone("manual");
    auto area_scaled = area.clone();
    area_scaled.scale(1/area_sum);

    horiz_contraction(manual,scl_x,area_scaled,comm);

    // Compare
    auto diff = diag_f.clone();
    diff.update(manual,-1,1);
    auto diff_norm = inf_norm(diff,&comm).as<Real>();
    REQUIRE_THAT (diff_norm, Catch::WithinAbs(0,tol));
  }

  SECTION ("vector_xz_masked") {
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "vec_xz");

    FieldLayout layout = grid->get_3d_vector_layout(LEV,2,"dim2");
    FieldIdentifier vec_xz_fid ("vec_xz" , layout, ekat::units::none, grid->name());
    Field vec_xz (vec_xz_fid , true);
    randomize_uniform(vec_xz,seed++);
    vec_xz.get_header().get_tracking().update_time_stamp(t0);
    auto& mask = vec_xz.create_valid_mask();
    randomize_uniform(mask,seed++);

    auto diag = diag_factory.create("HorizAvg", comm, params, grid);
    diag->set_input_field(vec_xz);
    diag->initialize();
    diag->compute_diagnostic(t0);
    auto diag_f = diag->get_diagnostic();
    REQUIRE (diag_f.has_valid_mask());

    auto& d_mask = diag_f.get_valid_mask();

    // Manual calculation: numerator = sum_col(area*f*mask), denominator = sum_col(area*mask)
    auto manual = diag_f.clone("manual");
    horiz_contraction(manual,vec_xz,area,comm);

    auto ones_xz = vec_xz.clone("ones");
    ones_xz.deep_copy(1);
    ones_xz.set_valid_mask(mask);
    auto manual_den = diag_f.clone("manual_den");
    horiz_contraction(manual_den,ones_xz,area,comm);

    // Expected valid_mask: 1 where denominator > 0
    auto tgt_mask = d_mask.clone("tgt_mask");
    compute_mask(manual_den,0,Comparison::NE,tgt_mask);
    REQUIRE (views_are_equal(tgt_mask,d_mask));

    // Compute manual average: numerator / denominator (where valid)
    manual.scale_inv(manual_den,tgt_mask);

    // Compare where valid; zero out fill_value entries before norm
    auto diff = diag_f.clone("diff");
    diff.update(manual,-1,1,d_mask);
    diff.deep_copy(0,d_mask,true);
    auto diff_norm = inf_norm(diff,&comm).as<Real>();
    REQUIRE_THAT (diff_norm, Catch::WithinAbs(0,tol));
  }

  SECTION ("vector_xz_fully_masked") {
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "vec_xz");

    FieldLayout layout = grid->get_3d_vector_layout(LEV,2,"dim2");
    FieldIdentifier vec_xz_fid ("vec_xz" , layout, ekat::units::none, grid->name());
    Field vec_xz (vec_xz_fid , true);
    randomize_uniform(vec_xz,seed++);
    vec_xz.get_header().get_tracking().update_time_stamp(t0);
    auto& mask = vec_xz.create_valid_mask();
    mask.deep_copy(0);

    auto diag = diag_factory.create("HorizAvg", comm, params, grid);
    diag->set_input_field(vec_xz);
    diag->initialize();
    diag->compute_diagnostic(t0);
    auto diag_f = diag->get_diagnostic();
    REQUIRE (diag_f.has_valid_mask());

    auto& d_mask = diag_f.get_valid_mask();

    auto tgt_mask = d_mask.clone("tgt_mask");
    mask.deep_copy(0);
    diag->compute_diagnostic(t0);
    tgt_mask.deep_copy(0);
    REQUIRE(views_are_equal(d_mask,tgt_mask));
  }

  SECTION ("constant_x") {
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "s_xz");

    FieldLayout layout = grid->get_3d_scalar_layout(LEV);
    FieldIdentifier scl_xz_fid ("s_xz" , layout, ekat::units::none, grid->name());
    Field scl_xz (scl_xz_fid , true);
    scl_xz.get_header().get_tracking().update_time_stamp(t0);

    Real c = 1.2345;
    scl_xz.deep_copy(c);

    auto diag = diag_factory.create("HorizAvg", comm, params, grid);
    diag->set_input_field(scl_xz);
    diag->initialize();
    diag->compute_diagnostic(t0);
    auto diag_f = diag->get_diagnostic();

    // Constant fields should have an avg equal to the constant(within rounding errors)
    auto manual = diag_f.clone("manual");
    manual.deep_copy(c);

    // Compare
    auto diff = diag_f.clone();
    diff.update(manual,-1,1);
    auto diff_norm = inf_norm(diff,&comm).as<Real>();
    REQUIRE_THAT (diff_norm, Catch::WithinAbs(0,tol));
  }
}

}  // namespace scream
