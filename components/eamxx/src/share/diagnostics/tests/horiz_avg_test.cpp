#include <catch2/catch.hpp>

#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

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
    diag->compute(t0);
    auto diag_f = diag->get();

    auto manual = diag_f.clone("manual");
    auto area_scaled = area.clone(CloneFlags::CopyData);
    area_scaled.scale(1/area_sum);

    horiz_contraction(manual,scl_x,area_scaled,comm);

    // Compare
    auto diff = diag_f.clone(CloneFlags::CopyData);
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
    diag->compute(t0);
    auto diag_f = diag->get();
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
    auto diff = diag_f.clone("diff", CloneFlags::CopyData);
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
    diag->compute(t0);
    auto diag_f = diag->get();
    REQUIRE (diag_f.has_valid_mask());

    auto& d_mask = diag_f.get_valid_mask();

    auto tgt_mask = d_mask.clone("tgt_mask");
    mask.deep_copy(0);
    diag->compute(t0);
    tgt_mask.deep_copy(0);
    REQUIRE(views_are_equal(d_mask,tgt_mask));
  }

  SECTION ("weights") {
    // As in xarray: weighting_method='none' is X.mean('col'), weight_name is
    // X.weighted(W).mean('col'), and weight_name='area' is the grid's cell area
    FieldLayout layout = grid->get_2d_scalar_layout();
    Field x (FieldIdentifier("x", layout, ekat::units::K, grid->name()), true);
    Field w (FieldIdentifier("w", layout, ekat::units::Pa, grid->name()), true);
    for (auto f : {&x,&w}) {
      f->get_header().get_tracking().update_time_stamp(t0);
    }
    randomize_uniform(x,seed++);
    randomize_uniform(w,seed++,0.5,2.0);

    auto make = [&](const std::string& key, const std::string& val) {
      ekat::ParameterList params;
      params.set("grid_name", grid->name());
      params.set<std::string>("field_name", "x");
      params.set<std::string>(key, val);
      return diag_factory.create("HorizAvg", comm, params, grid);
    };
    auto run = [&](const std::shared_ptr<AbstractDiagnostic>& d) {
      d->set_input_field(x);
      if (d->get_params().isParameter("weight_name") and
          d->get_params().get<std::string>("weight_name")=="w") {
        d->set_input_field(w);
      }
      d->initialize();
      d->compute(t0);
      return d->get();
    };
    auto check = [&](const Field& diag_f, const Field& weight) {
      auto num = diag_f.clone("num");
      auto den = diag_f.clone("den");
      auto ones = x.clone("ones");
      ones.deep_copy(1);
      horiz_contraction(num,x,weight,comm);
      horiz_contraction(den,ones,weight,comm);
      num.scale_inv(den);
      auto diff = diag_f.clone(CloneFlags::CopyData);
      diff.update(num,-1,1);
      REQUIRE_THAT (inf_norm(diff,&comm).as<Real>(), Catch::WithinAbs(0,tol));
    };

    // Unweighted: the plain average over all columns
    auto ones = w.clone("ones");
    ones.deep_copy(1);
    auto du = run(make("weighting_method","none"));
    REQUIRE (not du.has_valid_mask());
    check(du,ones);

    // 'area' is the grid's area, which is not an input, and the default
    auto da = make("weight_name","area");
    REQUIRE (da->get_input_fields_names()==std::list<std::string>{"x"});
    check(run(da),area);

    // A weight field: sum(w*x)/sum(w), masked only where the weights sum to 0
    auto dw = run(make("weight_name","w"));
    REQUIRE (dw.has_valid_mask());
    REQUIRE (dw.get_header().get_identifier().get_units()==ekat::units::K);
    dw.get_valid_mask().sync_to_host();
    REQUIRE (dw.get_valid_mask().get_internal_view_data<int,Host>()[0]==1);
    check(dw,w);

    w.deep_copy(0);
    auto dz = run(make("weight_name","w"));
    dz.sync_to_host();
    dz.get_valid_mask().sync_to_host();
    REQUIRE (dz.get_valid_mask().get_internal_view_data<int,Host>()[0]==0);
    REQUIRE (dz.get_internal_view_data<Real,Host>()[0]==constants::fill_value<Real>);

    // Rejected: a masked weight, a weight that is not (COL), both params, a bad method
    auto wm = w.clone("wm");
    wm.create_valid_mask();
    auto dm = make("weight_name","wm");
    dm->set_input_field(x);
    dm->set_input_field(wm);
    REQUIRE_THROWS (dm->initialize());
    Field w3 (FieldIdentifier("w3", grid->get_3d_scalar_layout(LEV), ekat::units::Pa, grid->name()), true);
    auto d3 = make("weight_name","w3");
    d3->set_input_field(x);
    d3->set_input_field(w3);
    REQUIRE_THROWS (d3->initialize());
    ekat::ParameterList both;
    both.set("grid_name", grid->name());
    both.set<std::string>("field_name", "x");
    both.set<std::string>("weight_name", "w");
    both.set<std::string>("weighting_method", "none");
    REQUIRE_THROWS (diag_factory.create("HorizAvg", comm, both, grid));
    REQUIRE_THROWS (make("weighting_method","dp"));
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
    diag->compute(t0);
    auto diag_f = diag->get();

    // Constant fields should have an avg equal to the constant(within rounding errors)
    auto manual = diag_f.clone("manual");
    manual.deep_copy(c);

    // Compare
    auto diff = diag_f.clone(CloneFlags::CopyData);
    diff.update(manual,-1,1);
    auto diff_norm = inf_norm(diff,&comm).as<Real>();
    REQUIRE_THAT (diff_norm, Catch::WithinAbs(0,tol));
  }
}

}  // namespace scream
