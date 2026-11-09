#include "catch2/catch.hpp"

#include "share/diagnostics/height_level_index.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_utils.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

namespace scream {

std::shared_ptr<HeightLevelIndex>
get_index_diag (const Field& z, std::shared_ptr<const AbstractGrid> grid,
                const std::string& layer, const Real height,
                const std::string& surf_ref)
{
  ekat::ParameterList params;
  params.set("surface_reference",surf_ref);
  params.set("height_value",std::to_string(height));
  params.set("height_units",std::string("m"));
  params.set("vertical_layer",layer);
  const auto& comm = grid->get_comm();
  auto diag = std::make_shared<HeightLevelIndex>(comm,params,grid);
  // The diag expects its input field to be named "z_<layer>" (sealevel) or
  // "height_<layer>" (surface); present the (single) input data under
  // whichever name is needed for this surf_ref, rather than requiring the
  // caller to keep separate, identically-valued fields around.
  const std::string z_name = std::string(surf_ref=="sealevel" ? "z" : "height") + "_" + layer;
  auto z_in = z.name()==z_name ? z : z.alias(z_name);
  diag->set_input_field(z_in);
  return diag;
}

TEST_CASE("height_level_index")
{
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);

  int ncols = 3;
  int nlevs = 10;
  auto grid = create_point_grid("physics",ncols,nlevs,comm);

  const Real z_top = 1000.0;
  const Real dz = z_top/nlevs;

  auto layout_mid = grid->get_3d_scalar_layout(LEV);
  FieldIdentifier z_mid_fid ("z_mid",layout_mid,ekat::units::m,grid->name());
  Field z_mid (z_mid_fid,true);

  auto z_mid_h = z_mid.get_view<Real**,Host>();
  for (int icol=0; icol<ncols; ++icol) {
    for (int ilev=0; ilev<nlevs; ++ilev) {
      // Descending with level index, as is standard for mid-level heights.
      z_mid_h(icol,ilev) = z_top - (ilev+0.5)*dz;
    }
  }
  util::TimeStamp t0 ({2022,1,1},{0,0,0});
  z_mid.get_header().get_tracking().update_time_stamp(t0);
  z_mid.sync_to_dev();

  // Compute the diag, and return its (host) index and range mask
  auto run = [&](const Real z_tgt, const std::string& surf_ref) {
    auto diag = get_index_diag(z_mid,grid,"mid",z_tgt,surf_ref);
    diag->initialize();
    diag->compute(t0);
    auto out = diag->get();
    out.sync_to_host();
    auto mask_f = out.get_valid_mask();
    mask_f.sync_to_host();
    return std::make_pair(out,mask_f);
  };

  // Check the same index and mask in all columns (the test heights are the same in all columns)
  auto check = [&](const Real z_tgt, const std::string& surf_ref,
                   const int k_expected, const int mask_expected) {
    auto [out,mask_f] = run(z_tgt,surf_ref);
    auto idx = out.get_view<const int*,Host>();
    auto mask = mask_f.get_view<const int*,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE(idx(icol)==k_expected);
      REQUIRE(mask(icol)==mask_expected);
    }
  };

  // The index is the first level with height below the target (z is descending
  // with level index), so the target is between the levels k-1 and k.
  auto check_in_range = [&](const std::string& surf_ref) {
    auto engine = scream::setup_random_test(&comm);
    std::uniform_real_distribution<Real> pdf(dz,z_top-dz);
    for (int itest=0; itest<10; ++itest) {
      Real z_tgt = pdf(engine);
      auto [out,mask_f] = run(z_tgt,surf_ref);
      auto idx = out.get_view<const int*,Host>();
      auto mask = mask_f.get_view<const int*,Host>();
      for (int icol=0; icol<ncols; ++icol) {
        const int k = idx(icol);
        REQUIRE(mask(icol)==1);
        REQUIRE(k>=1);
        REQUIRE(k<=nlevs-1);
        REQUIRE(z_mid_h(icol,k-1)>=z_tgt);
        REQUIRE(z_mid_h(icol,k)<z_tgt);
      }
    }
  };

  SECTION ("in_range_sealevel") {
    check_in_range("sealevel");
  }

  SECTION ("in_range_surface") {
    check_in_range("surface");
  }

  // A target above the top entry is always flagged as invalid, regardless
  // of the surface reference (matches PressureLevelIndex at the model top).
  SECTION ("top_invalid_surface") {
    check(2*z_top,"surface",0,0);
  }

  SECTION ("top_invalid_sealevel") {
    check(2*z_top,"sealevel",0,0);
  }

  // A target equal to the top entry is in range
  SECTION ("top_equal_is_valid") {
    check(z_mid_h(0,0),"sealevel",1,1);
    check(z_mid_h(0,0),"surface",1,1);
  }

  // A target below the bottom entry extrapolates for "above surface"
  // heights (always well defined near the surface): the index is nlevs, and the mask is 1...
  SECTION ("bottom_extrapolates_for_surface") {
    check(0,"surface",nlevs,1);
  }

  // ...but is flagged as invalid for "above sealevel" heights (e.g. a
  // target elevation below a mountain's surface is not well defined).
  SECTION ("bottom_invalid_for_sealevel") {
    check(0,"sealevel",nlevs,0);
  }
}

} // namespace scream
