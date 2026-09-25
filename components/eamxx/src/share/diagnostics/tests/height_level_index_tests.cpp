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

  auto check_in_range = [&](const std::string& surf_ref) {
    auto engine = scream::setup_random_test(&comm);
    std::uniform_real_distribution<Real> pdf(dz,z_top-dz);
    for (int itest=0; itest<10; ++itest) {
      Real z_tgt = pdf(engine);
      auto diag = get_index_diag(z_mid,grid,"mid",z_tgt,surf_ref);
      diag->initialize();
      diag->compute(t0);
      auto out = diag->get();
      out.sync_to_host();
      auto idx = out.get_view<const int**,Host>();
      for (int icol=0; icol<ncols; ++icol) {
        auto k0 = idx(icol,0);
        auto k1 = idx(icol,1);
        if (k0==k1) {
          REQUIRE(z_mid_h(icol,k0)==z_tgt);
        } else {
          REQUIRE(k1==k0+1);
          // z is descending with level index, so k0 (smaller index) has the larger z.
          REQUIRE(z_mid_h(icol,k0)>=z_tgt);
          REQUIRE(z_mid_h(icol,k1)<=z_tgt);
        }
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
    Real z_tgt = 2*z_top;
    auto diag = get_index_diag(z_mid,grid,"mid",z_tgt,"surface");
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

  SECTION ("top_invalid_sealevel") {
    Real z_tgt = 2*z_top;
    auto diag = get_index_diag(z_mid,grid,"mid",z_tgt,"sealevel");
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

  // A target below the bottom entry extrapolates for "above surface"
  // heights (always well defined near the surface)...
  SECTION ("bottom_extrapolates_for_surface") {
    Real z_tgt = 0;
    auto diag = get_index_diag(z_mid,grid,"mid",z_tgt,"surface");
    diag->initialize();
    diag->compute(t0);
    auto out = diag->get();
    out.sync_to_host();
    auto idx = out.get_view<const int**,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      REQUIRE(idx(icol,0)==nlevs-1);
      REQUIRE(idx(icol,1)==nlevs-1);
    }
  }

  // ...but is flagged as invalid for "above sealevel" heights (e.g. a
  // target elevation below a mountain's surface is not well defined).
  SECTION ("bottom_invalid_for_sealevel") {
    Real z_tgt = 0;
    auto diag = get_index_diag(z_mid,grid,"mid",z_tgt,"sealevel");
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

} // namespace scream
