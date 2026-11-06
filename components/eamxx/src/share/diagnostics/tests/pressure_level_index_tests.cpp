#include "catch2/catch.hpp"

#include "share/diagnostics/vertical_level_index.hpp"
#include "share/diagnostics/level_index_utils.hpp"
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

std::shared_ptr<VerticalLevelIndex>
get_index_diag (const Field& p, std::shared_ptr<const AbstractGrid> grid,
                const std::string& layer, const Real plevel)
{
  ekat::ParameterList params;
  params.set("vertical_coordinate",std::string("p"));
  params.set("value",std::to_string(plevel));
  params.set("units",std::string("Pa"));
  params.set("vertical_layer",layer);
  const auto& comm = grid->get_comm();
  auto diag = std::make_shared<VerticalLevelIndex>(comm,params,grid);
  diag->set_input_field(p);
  return diag;
}

TEST_CASE("pressure_level_in_Pa")
{
  // Different spellings of the same level must give the same value in Pa
  using namespace level_index_utils;
  const auto p1 = pressure_in_Pa("900","hPa");
  const auto p2 = pressure_in_Pa("900","mb");
  const auto p3 = pressure_in_Pa("90000","Pa");
  REQUIRE (p1==p3);
  REQUIRE (p2==p3);
  REQUIRE_THROWS (pressure_in_Pa("900","bar"));
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

  // Compute the diag, and return its (host) index and range mask
  auto run = [&](const Field& p, const Real plevel) {
    auto diag = get_index_diag(p,grid,"mid",plevel);
    diag->initialize();
    diag->compute(t0);
    auto out = diag->get();
    out.sync_to_host();
    auto mask_f = out.get_valid_mask();
    mask_f.sync_to_host();
    return std::make_pair(out,mask_f);
  };

  SECTION ("in_range") {
    // A target pressure strictly inside the column's range: the index is the first
    // level with pressure greater than the target, so the target is between k-1 and k.
    // Note: p_mid never actually reaches p_top/p_surf in any column (see
    // get_test_pres), so offset the bounds to guarantee we stay in range
    // for every column.
    auto engine = scream::setup_random_test(&comm);
    Real p_mid_bnds_dz = bnds.p_surf/nlevs;
    std::uniform_real_distribution<Real> pdf(bnds.p_top+p_mid_bnds_dz,bnds.p_surf-p_mid_bnds_dz);
    for (int itest=0; itest<10; ++itest) {
      Real plevel = std::round(pdf(engine));
      auto [out,mask_f] = run(p_mid,plevel);
      auto idx = out.get_view<const int*,Host>();
      auto mask = mask_f.get_view<const int*,Host>();
      for (int icol=0; icol<ncols; ++icol) {
        const int k = idx(icol);
        REQUIRE(mask(icol)==1);
        REQUIRE(k>=1);
        REQUIRE(k<=nlevs-1);
        REQUIRE(p_mid_h(icol,k-1)<=plevel);
        REQUIRE(p_mid_h(icol,k)>plevel);
      }
    }
  }

  SECTION ("out_of_range") {
    // Above the last level
    auto [out,mask_f] = run(p_mid,bnds.p_surf*2);
    auto idx = out.get_view<const int*,Host>();
    auto mask = mask_f.get_view<const int*,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE(mask(icol)==0);
      REQUIRE(idx(icol)==nlevs);
    }
  }

  SECTION ("edges") {
    // Use integer pressures, which survive the conversion to string: 100,200,...,1000 Pa
    FieldIdentifier p_fid ("p_mid",layout_mid,Pa,grid->name());
    Field p (p_fid,true);
    auto p_v = p.get_view<Real**,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      for (int ilev=0; ilev<nlevs; ++ilev) {
        p_v(icol,ilev) = 100.0*(ilev+1);
      }
    }
    p.get_header().get_tracking().update_time_stamp(t0);
    p.sync_to_dev();

    auto check = [&](const Real plevel, const int k_expected, const int mask_expected) {
      auto [out,mask_f] = run(p,plevel);
      auto idx = out.get_view<const int*,Host>();
      auto mask = mask_f.get_view<const int*,Host>();
      for (int icol=0; icol<ncols; ++icol) {
        REQUIRE(idx(icol)==k_expected);
        REQUIRE(mask(icol)==mask_expected);
      }
    };

    check(550,5,1);       // Between the levels 4 (500Pa) and 5 (600Pa)
    check(100,1,1);       // Equal to the first level: still in range
    check(1000,nlevs,1);  // Equal to the last level: in range, even if the index is nlevs
    check(50,0,0);        // Above the top
    check(1001,nlevs,0);  // Below the bottom (same index as the case above: the mask tells them apart)
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
