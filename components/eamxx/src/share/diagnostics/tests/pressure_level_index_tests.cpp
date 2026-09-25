#include "catch2/catch.hpp"

#include "share/diagnostics/pressure_level_index.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

struct PressureBnds
{
  Real p_top  = 10000.0;  //  100mb
  Real p_surf = 100000.0; // 1000mb
};

Real get_test_pres(const int col, const int lev, const int num_lev, const int num_cols);

std::shared_ptr<PressureLevelIndex>
get_index_diag (const Field& p, std::shared_ptr<const AbstractGrid> grid,
                const std::string& layer, const Real plevel)
{
  ekat::ParameterList params;
  params.set("pressure_value",std::to_string(plevel));
  params.set("pressure_units",std::string("Pa"));
  params.set("vertical_layer",layer);
  const auto& comm = grid->get_comm();
  auto diag = std::make_shared<PressureLevelIndex>(comm,params,grid);
  diag->set_input_field(p);
  return diag;
}

TEST_CASE("pressure_level_index")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  int ncols = 3;
  int nlevs = 10;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  using namespace ekat::units;
  auto layout_mid = grid->get_3d_scalar_layout(LEV);
  FieldIdentifier p_mid_fid ("p_mid",layout_mid,Pa,grid->name());
  Field p_mid (p_mid_fid,true);

  for (int icol=0; icol<ncols; ++icol) {
    for (int ilev=0; ilev<nlevs; ++ilev) {
      p_mid.get_view<Real**,Host>()(icol,ilev) = get_test_pres(icol,ilev,nlevs,ncols);
    }
  }
  util::TimeStamp t0 ({2022,1,1},{0,0,0});
  p_mid.get_header().get_tracking().update_time_stamp(t0);
  p_mid.sync_to_dev();

  auto p_mid_h = p_mid.get_view<const Real**,Host>();

  PressureBnds bnds;

  SECTION ("in_range") {
    // A target pressure strictly inside the column's range should give a
    // valid bracket, straddling the target pressure.
    // Note: p_mid never actually reaches p_top/p_surf in any column (see
    // get_test_pres), so offset the bounds to guarantee we stay in range
    // for every column.
    auto engine = scream::setup_random_test(&comm);
    Real p_mid_bnds_dz = bnds.p_surf/nlevs;
    std::uniform_real_distribution<Real> pdf(bnds.p_top+p_mid_bnds_dz,bnds.p_surf-p_mid_bnds_dz);
    for (int itest=0; itest<10; ++itest) {
      Real plevel = std::round(pdf(engine));
      auto diag = get_index_diag(p_mid,grid,"mid",plevel);
      diag->initialize();
      diag->compute(t0);
      auto out = diag->get();
      out.sync_to_host();
      auto idx = out.get_view<const int**,Host>();
      for (int icol=0; icol<ncols; ++icol) {
        auto k0 = idx(icol,0);
        auto k1 = idx(icol,1);
        // k0 and k1 are always distinct, valid, in-bounds indices.
        REQUIRE(k0>=0);
        REQUIRE(k1==k0+1);
        REQUIRE(p_mid_h(icol,k0)<=plevel);
        REQUIRE(p_mid_h(icol,k1)>=plevel);
      }
    }
  }

  SECTION ("out_of_range") {
    Real plevel = bnds.p_surf*2;
    auto diag = get_index_diag(p_mid,grid,"mid",plevel);
    diag->initialize();
    diag->compute(t0);
    auto out = diag->get();
    out.sync_to_host();
    auto idx = out.get_view<const int**,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE(idx(icol,0)==-1);
      REQUIRE(idx(icol,1)==-1);
    }
  }
}

Real get_test_pres(const int col, const int lev, const int num_lev, const int num_cols)
{
  PressureBnds bnds;
  Real p_surf = bnds.p_surf;
  Real p_top  = bnds.p_top;
  Real dp_dx  = p_top/(num_cols);
  Real dp_dz  = (p_surf-(p_top-(col+1)*dp_dx))/(num_lev);
  return (p_top - (col+1)*dp_dx) + lev*dp_dz;
}

} // namespace scream
