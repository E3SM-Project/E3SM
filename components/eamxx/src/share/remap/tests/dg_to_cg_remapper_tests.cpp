#include <catch2/catch.hpp>

#include "share/remap/dg_to_cg_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

constexpr int num_elems = 3;
constexpr int num_gp    = 4;
constexpr int num_levs  = 2*SCREAM_PACK_SIZE+1;
constexpr int vec_dim   = 3;

using gid_type = AbstractGrid::gid_type;

// The CG gid of the gp (ie,ip,jp) on this rank. Elements form a strip,
// where consecutive elements share a edge: ip=num_gp-1 on element ie is
// the same dof as ip=0 on element ie+1. Different ranks have disjoint dofs.
gid_type cg_gid (const int rank, const int ie, const int ip, const int jp)
{
  const int num_cg_per_rank = (num_elems*(num_gp-1)+1)*num_gp;
  return rank*num_cg_per_rank + (ie*(num_gp-1)+ip)*num_gp + jp;
}

// The value stored in a field at a given gid. The weights make sure that
// different components/levels contain different values
// (we assume icmp<10 and ilev<100)
Real value (const gid_type gid, const int icmp, const int ilev)
{
  // NOTE: these are integers (stored as Real), so that they are exact, no matter how the
  //       compiler contracts/vectorizes the arithmetic at different call sites
  return Real(gid)*1000 + icmp*100 + ilev;
}

std::shared_ptr<SEGrid> build_se_grid (const ekat::Comm& comm)
{
  auto se_grid = std::make_shared<SEGrid>("se",num_elems,num_gp,num_levs,comm);

  auto dg_gids   = se_grid->get_dofs_gids();
  auto cg_gids   = se_grid->get_cg_dofs_gids();
  auto elem_gids = se_grid->get_partitioned_dim_gids();
  auto lid2idx   = se_grid->get_lid_to_idx_map();

  auto h_dg    = dg_gids.get_view<gid_type*,Host>();
  auto h_cg    = cg_gids.get_view<gid_type*,Host>();
  auto h_elem  = elem_gids.get_view<gid_type*,Host>();
  auto h_l2i   = lid2idx.get_view<int**,Host>();

  const int rank = comm.rank();
  for (int ie=0; ie<num_elems; ++ie) {
    h_elem(ie) = ie + num_elems*rank;
    for (int ip=0; ip<num_gp; ++ip) {
      for (int jp=0; jp<num_gp; ++jp) {
        const int idof = ie*num_gp*num_gp + ip*num_gp + jp;
        h_dg(idof) = rank*num_elems*num_gp*num_gp + idof;
        h_cg(idof) = cg_gid(rank,ie,ip,jp);
        h_l2i(idof,0) = ie;
        h_l2i(idof,1) = ip;
        h_l2i(idof,2) = jp;
      }
    }
  }
  dg_gids.sync_to_dev();
  cg_gids.sync_to_dev();
  elem_gids.sync_to_dev();
  lid2idx.sync_to_dev();

  return se_grid;
}

// A point grid containing a subset of the unique dofs of the SE grid: the first row
// of each element, plus the last row of the strip (in a different order than the SE grid)
std::shared_ptr<PointGrid> build_pt_grid (const ekat::Comm& comm)
{
  std::vector<gid_type> gids;
  const int rank = comm.rank();
  for (int ie=num_elems-1; ie>=0; --ie) {
    for (int jp=0; jp<num_gp; ++jp) {
      gids.push_back(cg_gid(rank,ie,0,jp));
    }
  }
  // The last row of the last element is not the first row of any other element
  for (int jp=0; jp<num_gp; ++jp) {
    gids.push_back(cg_gid(rank,num_elems-1,num_gp-1,jp));
  }

  auto pt_grid = std::make_shared<PointGrid>("pt",gids.size(),num_levs,comm);
  auto pt_gids = pt_grid->get_dofs_gids();
  auto h_gids  = pt_gids.get_view<gid_type*,Host>();
  for (size_t i=0; i<gids.size(); ++i) {
    h_gids(i) = gids[i];
  }
  pt_gids.sync_to_dev();

  return pt_grid;
}

TEST_CASE ("dg_to_cg_remap")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto se_grid = build_se_grid(comm);
  auto pt_grid = build_pt_grid(comm);

  SECTION ("bad_input") {
    // Wrong grid types
    REQUIRE_THROWS (std::make_shared<DGToCGRemapper>(pt_grid,pt_grid));
    REQUIRE_THROWS (std::make_shared<DGToCGRemapper>(se_grid,se_grid));
  }

  SECTION ("remap") {
    const int ncols = pt_grid->get_num_local_dofs();
    const auto pt_gids = pt_grid->get_dofs_gids().get_view<const gid_type*,Host>();

    auto remap = std::make_shared<DGToCGRemapper>(se_grid,pt_grid);
    REQUIRE (remap->fwd_allowed());
    REQUIRE (not remap->bwd_allowed());

    auto make_se_field = [&](const std::string& name, const FieldLayout& layout) {
      Field f(FieldIdentifier(name,layout,none,se_grid->name()));
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      f.allocate_view();
      return f;
    };

    // SE fields, filled with values depending on the cg gid of each dof
    auto s2d = make_se_field("s2d",se_grid->get_2d_scalar_layout());
    auto v2d = make_se_field("v2d",se_grid->get_2d_vector_layout(vec_dim,"vec"));
    auto s3d = make_se_field("s3d",se_grid->get_3d_scalar_layout(LEV));
    auto v3d = make_se_field("v3d",se_grid->get_3d_vector_layout(LEV,vec_dim,"vec"));

    const auto se_gids = se_grid->get_cg_dofs_gids().get_view<const gid_type*,Host>();
    const auto l2i     = se_grid->get_lid_to_idx_map().get_view<const int**,Host>();
    {
      auto h_s2d = s2d.get_view<Real***,Host>();
      auto h_v2d = v2d.get_view<Real****,Host>();
      auto h_s3d = s3d.get_view<Real****,Host>();
      auto h_v3d = v3d.get_view<Real*****,Host>();
      for (int idof=0; idof<se_grid->get_num_local_dofs(); ++idof) {
        const int ie = l2i(idof,0);
        const int ip = l2i(idof,1);
        const int jp = l2i(idof,2);
        const auto gid = se_gids(idof);
        h_s2d(ie,ip,jp) = value(gid,0,0);
        for (int ic=0; ic<vec_dim; ++ic) {
          h_v2d(ie,ic,ip,jp) = value(gid,ic,0);
        }
        for (int il=0; il<num_levs; ++il) {
          h_s3d(ie,ip,jp,il) = value(gid,0,il);
          for (int ic=0; ic<vec_dim; ++ic) {
            h_v3d(ie,ic,ip,jp,il) = value(gid,ic,il);
          }
        }
      }
      for (auto f : {s2d,v2d,s3d,v3d}) {
        f.sync_to_dev();
      }
    }

    // The tgt fields are created by the remapper. The last field uses a component of
    // a SE vector field as src, to test that we correctly handle subfields.
    auto p_s2d = remap->register_field_from_src(s2d);
    auto p_v2d = remap->register_field_from_src(v2d);
    auto p_s3d = remap->register_field_from_src(s3d);
    auto p_v3d = remap->register_field_from_src(v3d);
    const int icmp = 1;
    auto p_v3d_cmp = remap->register_field_from_src(v3d.get_component(icmp));
    remap->registration_ends();

    // Poison the tgt fields
    for (auto f : {p_s2d,p_v2d,p_s3d,p_v3d,p_v3d_cmp}) {
      f.deep_copy(-1);
    }

    remap->remap_fwd();

    for (auto f : {p_s2d,p_v2d,p_s3d,p_v3d,p_v3d_cmp}) {
      f.sync_to_host();
    }

    auto h_p_s2d = p_s2d.get_view<const Real*,Host>();
    auto h_p_v2d = p_v2d.get_view<const Real**,Host>();
    auto h_p_s3d = p_s3d.get_view<const Real**,Host>();
    auto h_p_v3d = p_v3d.get_view<const Real***,Host>();
    auto h_p_cmp = p_v3d_cmp.get_view<const Real**,Host>();
    for (int icol=0; icol<ncols; ++icol) {
      const auto gid = pt_gids(icol);
      CHECK (h_p_s2d(icol)==value(gid,0,0));
      for (int ic=0; ic<vec_dim; ++ic) {
        CHECK (h_p_v2d(icol,ic)==value(gid,ic,0));
      }
      for (int il=0; il<num_levs; ++il) {
        CHECK (h_p_s3d(icol,il)==value(gid,0,il));
        CHECK (h_p_cmp(icol,il)==value(gid,icmp,il));
        for (int ic=0; ic<vec_dim; ++ic) {
          CHECK (h_p_v3d(icol,ic,il)==value(gid,ic,il));
        }
      }
    }

    // Backward remap is not supported
    REQUIRE_THROWS (remap->remap_bwd());
  }
}

} // namespace scream
