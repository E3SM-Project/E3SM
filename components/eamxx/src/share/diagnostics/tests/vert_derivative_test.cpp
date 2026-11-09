#include "catch2/catch.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/point_grid.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_universal_constants.hpp"

namespace scream {

// VertDerivative is xarray's differentiate (numpy's gradient, edge_order=1),
// applied in each column with the coordinate's own values.
TEST_CASE("vert_derivative") {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  util::TimeStamp t0({2024, 1, 1}, {0, 0, 0});

  auto &diag_factory = DiagnosticFactory::instance();
  register_diagnostics();

  const auto tol = std::numeric_limits<Real>::epsilon() * 1e4;

  auto make_diag = [&](const std::shared_ptr<const AbstractGrid>& grid,
                       const std::string& field, const std::string& coord) {
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    params.set<std::string>("field_name", field);
    params.set<std::string>("coordinate", coord);
    return diag_factory.create("VertDerivative", comm, params, grid);
  };

  SECTION("bad_diag") {
    auto grid = create_point_grid("physics",4*comm.size(),5,comm);
    ekat::ParameterList params;
    params.set("grid_name", grid->name());
    // No field_name
    REQUIRE_THROWS(diag_factory.create("VertDerivative", comm, params, grid));
    // No coordinate
    params.set<std::string>("field_name", "qc");
    REQUIRE_THROWS(diag_factory.create("VertDerivative", comm, params, grid));
    // With respect to itself
    params.set<std::string>("coordinate", "qc");
    REQUIRE_THROWS(diag_factory.create("VertDerivative", comm, params, grid));
  }

  SECTION("numpy_reference") {
    // Reference values from numpy/xarray:
    //   x = [1,2,4,7,11], f = [1,3,2,5,4]
    //   np.gradient(f,x)       == xr.DataArray(f,coords={'x':x}).differentiate('x')
    //                          == [2, 7/6, 1/10, 13/28, -1/4]
    //   np.gradient(f,x[::-1]) == [-1/2, -1/42, -23/30, 1/6, 1]
    // Column i uses the coordinate scaled by (i+1), so its derivative is the
    // reference divided by (i+1).
    const int nlevs = 5;
    auto grid = create_point_grid("physics",4*comm.size(),nlevs,comm);
    const int ncols = grid->get_num_local_dofs();
    const auto lay = grid->get_3d_scalar_layout(LEV);

    const std::vector<Real> x   = {1, 2, 4, 7, 11};
    const std::vector<Real> fv  = {1, 3, 2, 5, 4};
    const std::vector<Real> inc = {2, 7.0/6, 0.1, 13.0/28, -0.25};
    const std::vector<Real> dec = {-0.5, -1.0/42, -23.0/30, 1.0/6, 1};

    Field f (FieldIdentifier("qc",    lay, kg/kg, grid->name()));
    Field p (FieldIdentifier("p_mid", lay, Pa,    grid->name()));
    Field z (FieldIdentifier("z_mid", lay, m,     grid->name()));
    for (auto fld : {&f,&p,&z}) {
      fld->allocate_view();
      fld->get_header().get_tracking().update_time_stamp(t0);
    }
    auto fh = f.get_view<Real**,Host>();
    auto ph = p.get_view<Real**,Host>();
    auto zh = z.get_view<Real**,Host>();
    for (int i=0; i<ncols; ++i) {
      for (int k=0; k<nlevs; ++k) {
        fh(i,k) = fv[k];
        ph(i,k) = (i+1)*x[k];          // increasing downward, like p_mid
        zh(i,k) = (i+1)*x[nlevs-1-k];  // decreasing downward, like z_mid
      }
    }
    f.sync_to_dev(); p.sync_to_dev(); z.sync_to_dev();

    for (const auto& [coord,ref] : {std::pair{&p,&inc},std::pair{&z,&dec}}) {
      auto diag = make_diag(grid,"qc",coord->name());
      diag->set_input_field(f);
      diag->set_input_field(*coord);
      diag->initialize();
      diag->compute(t0);
      auto d = diag->get();
      REQUIRE (d.get_header().get_identifier().get_units()==kg/kg/coord->get_header().get_identifier().get_units());
      REQUIRE (not d.has_valid_mask());
      d.sync_to_host();
      auto dh = d.get_view<const Real**,Host>();
      for (int i=0; i<ncols; ++i) {
        for (int k=0; k<nlevs; ++k) {
          REQUIRE_THAT (dh(i,k), Catch::Matchers::WithinAbs((*ref)[k]/(i+1),tol));
        }
      }
    }
  }

  SECTION("exact_on_quadratics") {
    // On any spacing, numpy's interior stencil is exact for a quadratic, and the
    // one-sided end stencil gives the secant slope. With f = 3 + 2x + x^2/2:
    //   interior: f'(x_k) = 2 + x_k;   ends: 2 + (x_0+x_1)/2 and 2 + (x_{n-2}+x_{n-1})/2
    const int nlevs = 72;
    auto grid = create_point_grid("physics",8*comm.size(),nlevs,comm);
    const int ncols = grid->get_num_local_dofs();
    const auto lay = grid->get_3d_scalar_layout(LEV);

    Field f  (FieldIdentifier("T_mid", lay, K,  grid->name()));
    Field x  (FieldIdentifier("p_mid", lay, Pa, grid->name()));
    Field dx (FieldIdentifier("dx",    lay, Pa, grid->name()));
    for (auto fld : {&f,&x,&dx}) {
      fld->allocate_view();
      fld->get_header().get_tracking().update_time_stamp(t0);
    }
    int seed = get_random_test_seed();
    randomize_uniform(dx, seed++, 0.2, 2.0);
    dx.sync_to_host();
    auto dxh = dx.get_view<const Real**,Host>();
    auto xh  = x.get_view<Real**,Host>();
    auto fh  = f.get_view<Real**,Host>();
    for (int i=0; i<ncols; ++i) {
      Real xk = 1;
      for (int k=0; k<nlevs; ++k) {
        xk += dxh(i,k);
        xh(i,k) = xk;
        fh(i,k) = 3 + 2*xk + xk*xk/2;
      }
    }
    x.sync_to_dev(); f.sync_to_dev();

    auto expected = [&](const int i, const int k) -> Real {
      if (k==0)       return 2 + (xh(i,0)+xh(i,1))/2;
      if (k==nlevs-1) return 2 + (xh(i,nlevs-2)+xh(i,nlevs-1))/2;
      return 2 + xh(i,k);
    };

    auto diag = make_diag(grid,"T_mid","p_mid");
    diag->set_input_field(f);
    diag->set_input_field(x);
    diag->initialize();
    diag->compute(t0);
    auto d = diag->get();
    d.sync_to_host();
    auto dh = d.get_view<const Real**,Host>();
    for (int i=0; i<ncols; ++i) {
      for (int k=0; k<nlevs; ++k) {
        REQUIRE_THAT (dh(i,k), Catch::Matchers::WithinRel(expected(i,k),tol));
      }
    }

    // Masked input: a level is valid where every level its stencil reads is valid
    // (k-1,k,k+1 inside; k,k+1 at the top; k-1,k at the bottom), which is where
    // xarray would not propagate a NaN, and there the value is unchanged.
    f.create_valid_mask();
    randomize_uniform(f.get_valid_mask(),seed++);
    f.get_valid_mask().sync_to_host();
    auto fm = f.get_valid_mask().get_view<const int**,Host>();

    auto mdiag = make_diag(grid,"T_mid","p_mid");
    mdiag->set_input_field(f);
    mdiag->set_input_field(x);
    mdiag->initialize();
    mdiag->compute(t0);
    auto md = mdiag->get();
    REQUIRE (md.has_valid_mask());
    REQUIRE (md.get_header().may_be_filled());
    md.sync_to_host();
    md.get_valid_mask().sync_to_host();
    auto mdh = md.get_view<const Real**,Host>();
    auto mdm = md.get_valid_mask().get_view<const int**,Host>();
    int nvalid = 0, ninvalid = 0;
    for (int i=0; i<ncols; ++i) {
      for (int k=0; k<nlevs; ++k) {
        const int km = k==0 ? 0 : k-1;
        const int kp = k==nlevs-1 ? nlevs-1 : k+1;
        int valid = 1;
        for (int j=km; j<=kp; ++j) valid *= fm(i,j)!=0 ? 1 : 0;
        REQUIRE ((mdm(i,k)!=0)==(valid==1));
        if (valid) {
          ++nvalid;
          REQUIRE (mdh(i,k)==dh(i,k));
        } else {
          ++ninvalid;
          REQUIRE (mdh(i,k)==constants::fill_value<Real>);
        }
      }
    }
    // Make sure the mask test could have failed both ways
    REQUIRE (nvalid>0);
    REQUIRE (ninvalid>0);
  }

  SECTION("layout_checks") {
    auto grid = create_point_grid("physics",4*comm.size(),5,comm);
    Field f (FieldIdentifier("qc",    grid->get_3d_scalar_layout(LEV), kg/kg, grid->name()));
    Field s (FieldIdentifier("ps",    grid->get_2d_scalar_layout(),    Pa,    grid->name()));
    Field pi(FieldIdentifier("p_int", grid->get_3d_scalar_layout(ILEV),Pa,    grid->name()));
    for (auto fld : {&f,&s,&pi}) {
      fld->allocate_view();
      fld->get_header().get_tracking().update_time_stamp(t0);
    }
    // Coordinate on a different layout
    auto d1 = make_diag(grid,"qc","p_int");
    d1->set_input_field(f);
    d1->set_input_field(pi);
    REQUIRE_THROWS (d1->initialize());
    // Field without a level dimension
    auto d2 = make_diag(grid,"ps","qc");
    d2->set_input_field(s);
    d2->set_input_field(f);
    REQUIRE_THROWS (d2->initialize());
  }
}

} // namespace scream
