#include <catch2/catch.hpp>

#include <cmath>
#include <limits>

#include "share/util/eamxx_universal_constants.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {


TEST_CASE("vert_contract") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  // A world comm
  ekat::Comm comm(MPI_COMM_WORLD);

  // A time stamp
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  // Create a grids manager
  const int nlevs = 30;
  const int dim5  = 5;
  const int ncols = 95;

  auto grid = create_point_grid("physics",ncols*comm.size(),nlevs,comm);

  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  int seed = get_random_test_seed();

  // Shared field layouts
  FieldLayout scalar2d_layout = grid->get_3d_scalar_layout(LEV);
  FieldLayout scalar3d_layout = grid->get_3d_vector_layout(LEV,dim5,"dim5");

  SECTION("bad_input") {
    ekat::ParameterList params;
    // No field_name -> throws in constructor
    REQUIRE_THROWS(diag_factory.create("VertContract", comm, params, grid));

    params.set<std::string>("field_name", "qc");
    params.set("grid_name", grid->name());
    // No contract_method -> throws in constructor
    REQUIRE_THROWS(diag_factory.create("VertContract", comm, params, grid));

    params.set<std::string>("contract_method", "xyz");
    // Invalid contract_method -> throws in constructor
    REQUIRE_THROWS(diag_factory.create("VertContract", comm, params, grid));
  }

  SECTION("dp_weighted_sum") {
    // Create and randomize a 3D (COL x CMP x LEV) input field
    Field fin(FieldIdentifier("qc", scalar3d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    Field dp(FieldIdentifier("pseudo_density", scalar2d_layout, Pa, grid->name()));
    dp.allocate_view();
    dp.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(dp, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "sum");
    params.set<std::string>("weighting_method", "dp");
    auto diag = diag_factory.create("VertContract", comm, params, grid);
    diag->set_input_field(fin);
    diag->set_input_field(dp);
    diag->initialize();
    diag->compute(t0);
    auto diag_f = diag->get();

    // Manual reference: ref(col,cmp) = sum_lev( (dp/g)(col,lev) * f(col,cmp,lev) )
    auto dp_scaled = dp.clone("dp_scaled", CloneFlags::CopyData);
    dp_scaled.scale(1 / physics::Constants<Real>::gravit.value);
    auto ref = diag_f.clone("ref");
    vert_contraction(ref, fin, dp_scaled);

    REQUIRE(views_are_equal(diag_f, ref));
  }

  SECTION("dp_weighted_avg") {
    // Create and randomize a 2D (COL x LEV) input field
    Field fin(FieldIdentifier("qc", scalar2d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    Field dp(FieldIdentifier("pseudo_density", scalar2d_layout, Pa, grid->name()));
    dp.allocate_view();
    dp.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(dp, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "avg");
    params.set<std::string>("weighting_method", "dp");
    auto diag = diag_factory.create("VertContract", comm, params, grid);
    diag->set_input_field(fin);
    diag->set_input_field(dp);
    diag->initialize();
    diag->compute(t0);
    auto diag_f = diag->get();

    // Manual reference:
    //   num(col) = sum_lev( (dp/g)(col,lev) * f(col,lev) )
    //   den(col) = sum_lev( (dp/g)(col,lev) )
    //   ref(col) = num(col) / den(col)
    auto dp_scaled = dp.clone("dp_scaled", CloneFlags::CopyData);
    dp_scaled.scale(1 / physics::Constants<Real>::gravit.value);
    auto dp_ones = dp_scaled.clone("ones");
    dp_ones.deep_copy(1);
    auto ref = diag_f.clone("ref");
    auto den = diag_f.clone("den");
    vert_contraction(ref, fin, dp_scaled);
    vert_contraction(den, dp_scaled, dp_ones);
    ref.scale_inv(den);

    REQUIRE(views_are_equal(diag_f, ref));
  }

  SECTION("unweighted_sum") {
    // Create and randomize a 2D (COL x LEV) input field
    Field fin(FieldIdentifier("qc", scalar2d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "sum");
    params.set<std::string>("weighting_method", "none");
    auto diag = diag_factory.create("VertContract", comm, params, grid);
    diag->set_input_field(fin);
    diag->initialize();
    diag->compute(t0);
    auto diag_f = diag->get();

    // Manual reference: ref(col) = sum_lev( f(col,lev) )
    // Use a 1D all-ones weight (same as the diagnostic)
    FieldLayout ones_lev_layout{{LEV}, {nlevs}};
    Field ones_lev(FieldIdentifier("ones_lev", ones_lev_layout, none, grid->name()), true);
    ones_lev.deep_copy(1);
    auto ref = diag_f.clone("ref");
    vert_contraction(ref, fin, ones_lev);

    REQUIRE(views_are_equal(diag_f, ref));
  }

  SECTION("masked_sum") {
    // Column 0 has no valid level, column 1 has its first three levels invalid
    Field fin(FieldIdentifier("qc", scalar2d_layout, kg/kg, grid->name()), true);
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);
    auto& fmask = fin.create_valid_mask(Field::MaskInit::Valid);
    fmask.sync_to_host();
    auto fm_h = fmask.get_view<int**,Host>();
    for (int k=0; k<nlevs; ++k) fm_h(0,k) = 0;
    for (int k=0; k<3; ++k)     fm_h(1,k) = 0;
    fmask.sync_to_dev();

    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "sum");
    params.set<std::string>("weighting_method", "none");
    auto diag = diag_factory.create("VertContract", comm, params, grid);
    diag->set_input_field(fin);
    diag->initialize();
    diag->compute(t0);
    auto diag_f = diag->get();

    // The sum of a masked field is masked where no level is valid
    REQUIRE(diag_f.has_valid_mask());

    diag_f.sync_to_host();
    diag_f.get_valid_mask().sync_to_host();
    fin.sync_to_host();
    auto d_h = diag_f.get_view<const Real*,Host>();
    auto m_h = diag_f.get_valid_mask().get_view<const int*,Host>();
    auto f_h = fin.get_view<const Real**,Host>();
    const int nlocal = grid->get_num_local_dofs();
    for (int i=0; i<nlocal; ++i) {
      if (i==0) {
        REQUIRE(m_h(i)==0);
      } else {
        Real ref = 0;
        for (int k=(i==1 ? 3 : 0); k<nlevs; ++k) ref += f_h(i,k);
        REQUIRE(m_h(i)==1);
        REQUIRE(d_h(i)==Approx(ref).epsilon(100*std::numeric_limits<Real>::epsilon()));
      }
    }
  }

  SECTION("weighted_by_field") {
    // X.sum/mean('lev',weights=W), against xarray (2025.10.1) on the same data,
    // with NaN where X is masked:
    //   x = [[1,nan,3],[nan,nan,nan],[2,4,nan],[5,6,7]]
    //   w = [[1,2,3],  [1,1,1],      [0,1,0],  [0,0,0]]
    //   x.weighted(w).sum('lev')  == [10, 0, 4, 0]
    //   x.weighted(w).mean('lev') == [2.5, nan, 4, nan]
    // Every rank holds the same four columns.
    constexpr auto fv = constants::fill_value<Real>;
    const int nc = 4, nl = 3;
    auto g4 = create_point_grid("physics",nc*comm.size(),nl,comm);
    const auto lay = g4->get_3d_scalar_layout(LEV);

    const Real xv[nc][nl] = {{1,fv,3},{fv,fv,fv},{2,4,fv},{5,6,7}};
    const int  xm[nc][nl] = {{1,0,1},{0,0,0},{1,1,0},{1,1,1}};
    const Real wv[nc][nl] = {{1,2,3},{1,1,1},{0,1,0},{0,0,0}};

    Field x(FieldIdentifier("x", lay, kg/kg, g4->name()));
    Field w(FieldIdentifier("w", lay, Pa,    g4->name()));
    x.allocate_view();
    w.allocate_view();
    x.create_valid_mask();
    auto xh = x.get_view<Real**,Host>();
    auto mh = x.get_valid_mask().get_view<int**,Host>();
    auto wh = w.get_view<Real**,Host>();
    for (int i=0; i<nc; ++i) {
      for (int k=0; k<nl; ++k) {
        xh(i,k) = xv[i][k];
        mh(i,k) = xm[i][k];
        wh(i,k) = wv[i][k];
      }
    }
    x.sync_to_dev(); x.get_valid_mask().sync_to_dev(); w.sync_to_dev();
    for (auto fld : {&x,&w}) {
      fld->get_header().get_tracking().update_time_stamp(t0);
    }

    auto make = [&](const std::string& method, const std::string& field, const std::string& wname) {
      ekat::ParameterList params;
      params.set("grid_name", g4->name());
      params.set<std::string>("field_name", field);
      params.set<std::string>("contract_method", method);
      params.set<std::string>("weight_name", wname);
      return diag_factory.create("VertContract", comm, params, g4);
    };

    // sum: invalid entries skipped, all-invalid column -> 0, no mask
    auto ds = make("sum","x","w");
    ds->set_input_field(x);
    ds->set_input_field(w);
    ds->initialize();
    ds->compute(t0);
    auto s = ds->get();
    REQUIRE (not s.has_valid_mask());
    REQUIRE (s.get_header().get_identifier().get_units()==kg/kg*Pa);
    s.sync_to_host();
    auto sh = s.get_view<const Real*,Host>();
    const Real sref[nc] = {10, 0, 4, 0};
    for (int i=0; i<nc; ++i) {
      REQUIRE (sh(i)==sref[i]);
    }

    // mean: divided by the weights of the valid entries; invalid where those sum to 0,
    // whether because no entry is valid (col 1) or because the weights vanish (col 3)
    auto dm = make("avg","x","w");
    dm->set_input_field(x);
    dm->set_input_field(w);
    dm->initialize();
    dm->compute(t0);
    auto m = dm->get();
    REQUIRE (m.has_valid_mask());
    REQUIRE (m.get_header().get_identifier().get_units()==kg/kg);
    m.sync_to_host();
    m.get_valid_mask().sync_to_host();
    auto mvh = m.get_view<const Real*,Host>();
    auto mmh = m.get_valid_mask().get_view<const int*,Host>();
    const int  mmask[nc] = {1, 0, 1, 0};
    const Real mref[nc]  = {2.5, 0, 4, 0};
    for (int i=0; i<nc; ++i) {
      REQUIRE (mmh(i)==mmask[i]);
      if (mmask[i]==1) {
        REQUIRE (mvh(i)==mref[i]);
      }
    }

    // Unmasked field: the mean is still invalid where the weights vanish (col 3)
    Field y(FieldIdentifier("y", lay, kg/kg, g4->name()));
    y.allocate_view();
    y.deep_copy(2);
    y.get_header().get_tracking().update_time_stamp(t0);
    auto dy = make("avg","y","w");
    dy->set_input_field(y);
    dy->set_input_field(w);
    dy->initialize();
    dy->compute(t0);
    auto my = dy->get();
    REQUIRE (my.has_valid_mask());
    my.sync_to_host();
    my.get_valid_mask().sync_to_host();
    auto myh = my.get_view<const Real*,Host>();
    auto mym = my.get_valid_mask().get_view<const int*,Host>();
    for (int i=0; i<nc; ++i) {
      REQUIRE (mym(i)==(i==3 ? 0 : 1));
      if (i!=3) {
        REQUIRE (myh(i)==Real(2));
      }
    }

    // Levels with zero weight are not read: X may be NaN there (xarray: NaN is
    // missing), and the sums are those of the finite values
    constexpr auto nan = std::numeric_limits<Real>::quiet_NaN();
    Field z(FieldIdentifier("z", lay, kg/kg, g4->name()));
    z.allocate_view();
    auto zh = z.get_view<Real**,Host>();
    const Real zv[nc][nl] = {{1,2,3},{1,1,1},{nan,4,nan},{nan,nan,nan}};
    for (int i=0; i<nc; ++i) {
      for (int k=0; k<nl; ++k) {
        zh(i,k) = zv[i][k];
      }
    }
    z.sync_to_dev();
    z.get_header().get_tracking().update_time_stamp(t0);
    const Real zsum[nc]  = {14, 3, 4, 0};
    const Real zmean[nc] = {14.0/6, 1, 4, fv};
    for (const std::string method : {"sum","avg"}) {
      auto dz = make(method,"z","w");
      dz->set_input_field(z);
      dz->set_input_field(w);
      dz->initialize();
      dz->compute(t0);
      auto o = dz->get();
      o.sync_to_host();
      auto oh = o.get_view<const Real*,Host>();
      for (int i=0; i<nc; ++i) {
        const Real ref = method=="sum" ? zsum[i] : zmean[i];
        REQUIRE (std::isfinite(oh(i)));
        REQUIRE (oh(i)==Approx(ref));
      }
    }

    // A field weighted by itself is one input: w.mean('lev',weights=w) = sum(w^2)/sum(w)
    auto dww = make("avg","w","w");
    dww->set_input_field(w);
    dww->initialize();
    dww->compute(t0);
    auto mww = dww->get();
    mww.sync_to_host();
    REQUIRE (mww.get_view<const Real*,Host>()(0)==Approx(14.0/6));

    // A masked weight is rejected (xarray: weights cannot be missing)
    auto dbad = make("sum","y","x");
    dbad->set_input_field(y);
    dbad->set_input_field(x);
    REQUIRE_THROWS (dbad->initialize());

    // A weight with components is rejected: (lev) or the field's layout without CMP
    Field wc(FieldIdentifier("wc", g4->get_3d_vector_layout(LEV,2,"dim2"), Pa, g4->name()));
    wc.allocate_view();
    wc.get_header().get_tracking().update_time_stamp(t0);
    auto dwc = make("sum","y","wc");
    dwc->set_input_field(y);
    dwc->set_input_field(wc);
    REQUIRE_THROWS (dwc->initialize());

    // weight_name replaces weighting_method; not both
    ekat::ParameterList both;
    both.set("grid_name", g4->name());
    both.set<std::string>("field_name", "y");
    both.set<std::string>("contract_method", "sum");
    both.set<std::string>("weighting_method", "dp");
    both.set<std::string>("weight_name", "w");
    REQUIRE_THROWS (diag_factory.create("VertContract", comm, both, g4));
  }

  SECTION("unweighted_avg") {
    // Create and randomize a 3D (COL x CMP x LEV) input field
    Field fin(FieldIdentifier("qc", scalar3d_layout, kg/kg, grid->name()));
    fin.allocate_view();
    fin.get_header().get_tracking().update_time_stamp(t0);
    randomize_uniform(fin, seed++, 1, 200);

    // Create and set up the diagnostic
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", "qc");
    params.set<std::string>("contract_method", "avg");
    params.set<std::string>("weighting_method", "none");
    auto diag = diag_factory.create("VertContract", comm, params, grid);
    diag->set_input_field(fin);
    diag->initialize();
    diag->compute(t0);
    auto diag_f = diag->get();

    // Manual reference:
    //   num(col,cmp) = sum_lev( f(col,cmp,lev) )
    //   den(col,cmp) = sum_lev( 1 ) = nlevs
    //   ref(col,cmp) = num(col,cmp) / den(col,cmp)
    // Use a 1D all-ones weight (same as the diagnostic)
    FieldLayout ones_lev_layout{{LEV}, {nlevs}};
    Field ones_lev(FieldIdentifier("ones_lev", ones_lev_layout,
                                   none, grid->name()), true);
    ones_lev.deep_copy(1);
    auto fin_ones = fin.clone("fin_ones");
    fin_ones.deep_copy(1);
    auto ref = diag_f.clone("ref");
    auto den = diag_f.clone("den");
    vert_contraction(ref, fin, ones_lev);
    vert_contraction(den, fin_ones, ones_lev);
    ref.scale_inv(den);

    REQUIRE(views_are_equal(diag_f, ref));
  }
}

}  // namespace scream
