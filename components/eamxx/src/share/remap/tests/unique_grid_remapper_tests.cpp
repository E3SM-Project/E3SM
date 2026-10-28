#include <catch2/catch.hpp>

#include "share/remap/unique_grid_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

constexpr int num_elems = 3;
constexpr int num_gp    = 4;
constexpr int num_levs  = 2*SCREAM_PACK_SIZE+1;
constexpr int vec_dim   = 3;

using gid_type = AbstractGrid::gid_type;

// The unique gid of the gp (ie,ip,jp) on this rank. Elements form a strip (also across ranks),
// where consecutive elements share an edge: ip=num_gp-1 on element ie is the same dof as ip=0
// on element ie+1. The last element of rank r shares an edge with the first element of rank r+1.
// A "row" is a set of num_gp dofs (jp=0..num_gp-1) with the same ip.
constexpr int rows_per_rank = num_elems*(num_gp-1);
gid_type unique_gid (const int rank, const int ie, const int ip, const int jp)
{
  const int row = rank*rows_per_rank + ie*(num_gp-1) + ip;
  return row*num_gp + jp;
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

std::shared_ptr<SEGrid> build_se_grid (const ekat::Comm& comm, Field& unique_gids)
{
  using namespace ShortFieldTagsNames;
  auto se_grid = std::make_shared<SEGrid>("se",num_elems,num_gp,num_levs,comm);

  auto dofs_gids   = se_grid->get_dofs_gids();
  unique_gids = Field(FieldIdentifier("unique_grid_gids",FieldLayout({CMP},{se_grid->get_num_local_dofs()}),
                                      ekat::units::none,se_grid->name(),DataType::IntType));
  unique_gids.allocate_view();
  auto elem_gids = se_grid->get_partitioned_dim_gids();
  auto lid2idx   = se_grid->get_lid_to_idx_map();

  auto h_dofs    = dofs_gids.get_view<gid_type*,Host>();
  auto h_unique    = unique_gids.get_view<gid_type*,Host>();
  auto h_elem  = elem_gids.get_view<gid_type*,Host>();
  auto h_l2i   = lid2idx.get_view<int**,Host>();

  const int rank = comm.rank();
  for (int ie=0; ie<num_elems; ++ie) {
    h_elem(ie) = ie + num_elems*rank;
    for (int ip=0; ip<num_gp; ++ip) {
      for (int jp=0; jp<num_gp; ++jp) {
        const int idof = ie*num_gp*num_gp + ip*num_gp + jp;
        h_dofs(idof) = rank*num_elems*num_gp*num_gp + idof;
        h_unique(idof) = unique_gid(rank,ie,ip,jp);
        h_l2i(idof,0) = ie;
        h_l2i(idof,1) = ip;
        h_l2i(idof,2) = jp;
      }
    }
  }
  dofs_gids.sync_to_dev();
  unique_gids.sync_to_dev();
  elem_gids.sync_to_dev();
  lid2idx.sync_to_dev();

  return se_grid;
}

// A point grid containing the unique dofs of the SE grid owned by this rank (in a different order than the SE grid).
// Each rank owns the rows from its first one to the one before the first row of the next rank
// (the last rank owns the last row too). Hence, each rank needs, for the CG->DG remap, the row owned by the next rank.
// If subset=true, only half of the rows are kept (the grid is then not enough for CG->DG).
std::shared_ptr<PointGrid> build_pt_grid (const ekat::Comm& comm, const bool subset = false)
{
  std::vector<gid_type> gids;
  const int rank = comm.rank();
  const bool last = rank==comm.size()-1;
  const int nrows = rows_per_rank + (last ? 1 : 0);
  for (int irow=nrows-1; irow>=0; --irow) {
    if (subset and irow%2==1) {
      continue;
    }
    for (int jp=0; jp<num_gp; ++jp) {
      gids.push_back((rank*rows_per_rank+irow)*num_gp + jp);
    }
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

TEST_CASE ("unique_grid_remap")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);

  Field se_unique_gids;
  auto se_grid = build_se_grid(comm,se_unique_gids);
  auto pt_grid = build_pt_grid(comm);
  auto pt_subset_grid = build_pt_grid(comm,true);

  SECTION ("bad_input") {
    // The src grid is already unique
    REQUIRE_THROWS (std::make_shared<UniqueGridRemapper>(pt_grid));
    REQUIRE_THROWS (std::make_shared<UniqueGridRemapper>(se_grid));
  }

  // NOTE: we use CHECK (rather than REQUIRE) for the data checks, so that all ranks keep going in case of failure
  //       on a subset of ranks. Otherwise, Catch would re-run the test case on the failing ranks only, and since
  //       the grids creation is collective, the test would hang.
  SECTION ("dg_to_cg") {
    const int ncols = pt_subset_grid->get_num_local_dofs();
    const auto pt_gids = pt_subset_grid->get_dofs_gids().get_view<const gid_type*,Host>();

    se_grid->set_unique_grid(pt_subset_grid,se_unique_gids);
    auto remap = std::make_shared<UniqueGridRemapper>(se_grid);
    REQUIRE (remap->fwd_allowed());
    REQUIRE (remap->bwd_allowed());

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

    const auto se_gids = se_grid->get_unique_grid_gids().get_view<const gid_type*,Host>();
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
  }

  SECTION ("cg_to_dg") {
    // Note: this requires comm among ranks
    se_grid->set_unique_grid(pt_grid,se_unique_gids);
    auto remap = std::make_shared<UniqueGridRemapper>(se_grid);

    auto make_se_field = [&](const std::string& name, const FieldLayout& layout) {
      Field f(FieldIdentifier(name,layout,none,se_grid->name()));
      f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
      f.allocate_view();
      f.deep_copy(-1);
      return f;
    };

    auto s2d = make_se_field("s2d",se_grid->get_2d_scalar_layout());
    auto v2d = make_se_field("v2d",se_grid->get_2d_vector_layout(vec_dim,"vec"));
    auto s3d = make_se_field("s3d",se_grid->get_3d_scalar_layout(LEV));
    auto v3d = make_se_field("v3d",se_grid->get_3d_vector_layout(LEV,vec_dim,"vec"));
    // A SE field with a component that is also remapped, to test that we handle subfields
    auto v3d_2 = make_se_field("v3d_2",se_grid->get_3d_vector_layout(LEV,vec_dim,"vec"));
    const int icmp = 1;

    // The CG fields are created by the remapper
    auto p_s2d = remap->register_field_from_src(s2d);
    auto p_v2d = remap->register_field_from_src(v2d);
    auto p_s3d = remap->register_field_from_src(s3d);
    auto p_v3d = remap->register_field_from_src(v3d);
    auto p_cmp = remap->register_field_from_src(v3d_2.get_component(icmp));
    remap->registration_ends();

    // Fill the CG fields with values depending on the gid
    const int ncols = pt_grid->get_num_local_dofs();
    const auto pt_gids = pt_grid->get_dofs_gids().get_view<const gid_type*,Host>();
    {
      auto h_s2d = p_s2d.get_view<Real*,Host>();
      auto h_v2d = p_v2d.get_view<Real**,Host>();
      auto h_s3d = p_s3d.get_view<Real**,Host>();
      auto h_v3d = p_v3d.get_view<Real***,Host>();
      auto h_cmp = p_cmp.get_view<Real**,Host>();
      for (int icol=0; icol<ncols; ++icol) {
        const auto gid = pt_gids(icol);
        h_s2d(icol) = value(gid,0,0);
        for (int ic=0; ic<vec_dim; ++ic) {
          h_v2d(icol,ic) = value(gid,ic,0);
        }
        for (int il=0; il<num_levs; ++il) {
          h_s3d(icol,il) = value(gid,0,il);
          h_cmp(icol,il) = value(gid,icmp,il);
          for (int ic=0; ic<vec_dim; ++ic) {
            h_v3d(icol,ic,il) = value(gid,ic,il);
          }
        }
      }
      for (auto f : {p_s2d,p_v2d,p_s3d,p_v3d,p_cmp}) {
        f.sync_to_dev();
      }
    }

    remap->remap_bwd();

    for (auto f : {s2d,v2d,s3d,v3d,v3d_2}) {
      f.sync_to_host();
    }

    // All dofs (including the ones on the edges, and the ones whose CG column is owned by another rank) must be set
    const auto se_gids = se_grid->get_unique_grid_gids().get_view<const gid_type*,Host>();
    const auto l2i     = se_grid->get_lid_to_idx_map().get_view<const int**,Host>();
    auto h_s2d = s2d.get_view<const Real***,Host>();
    auto h_v2d = v2d.get_view<const Real****,Host>();
    auto h_s3d = s3d.get_view<const Real****,Host>();
    auto h_v3d = v3d.get_view<const Real*****,Host>();
    auto h_v3d_2 = v3d_2.get_view<const Real*****,Host>();
    for (int idof=0; idof<se_grid->get_num_local_dofs(); ++idof) {
      const int ie = l2i(idof,0);
      const int ip = l2i(idof,1);
      const int jp = l2i(idof,2);
      const auto gid = se_gids(idof);
      CHECK (h_s2d(ie,ip,jp)==value(gid,0,0));
      for (int ic=0; ic<vec_dim; ++ic) {
        CHECK (h_v2d(ie,ic,ip,jp)==value(gid,ic,0));
      }
      for (int il=0; il<num_levs; ++il) {
        CHECK (h_s3d(ie,ip,jp,il)==value(gid,0,il));
        for (int ic=0; ic<vec_dim; ++ic) {
          CHECK (h_v3d(ie,ic,ip,jp,il)==value(gid,ic,il));
          // Only one component of v3d_2 was remapped. The others must be untouched
          CHECK (h_v3d_2(ie,ic,ip,jp,il)==(ic==icmp ? value(gid,icmp,il) : -1));
        }
      }
    }
  }
}

TEST_CASE ("unique_grid_remap_1d_src")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  const int rank = comm.rank();

  // Per rank, n unique dofs. The strip has n+1 dofs, with the last one a replica of the first
  const int n = 5;

  auto strip = std::make_shared<PointGrid>("strip",n+1,num_levs,comm);
  Field strip_unique_gids (FieldIdentifier("unique_grid_gids",FieldLayout({CMP},{n+1}),none,strip->name(),DataType::IntType));
  strip_unique_gids.allocate_view();
  auto unique = std::make_shared<PointGrid>("unique",n,num_levs,comm);
  {
    auto dofs = strip->get_dofs_gids();
    auto map  = strip_unique_gids;
    auto h_dofs = dofs.get_view<gid_type*,Host>();
    auto h_map  = map.get_view<gid_type*,Host>();
    for (int i=0; i<=n; ++i) {
      h_dofs(i) = rank*(n+1) + i;
      h_map(i)  = rank*n + (i%n);
    }
    dofs.sync_to_dev();
    map.sync_to_dev();

    // Unique dofs, in reverse order
    auto udofs = unique->get_dofs_gids();
    auto h_udofs = udofs.get_view<gid_type*,Host>();
    for (int i=0; i<n; ++i) {
      h_udofs(i) = rank*n + (n-1-i);
    }
    udofs.sync_to_dev();
  }
  strip->set_unique_grid(unique,strip_unique_gids);

  auto remap = std::make_shared<UniqueGridRemapper>(strip);
  REQUIRE (remap->fwd_allowed());

  auto make_field = [&](const std::string& name, const FieldLayout& layout) {
    Field f(FieldIdentifier(name,layout,none,strip->name()));
    f.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
    f.allocate_view();
    return f;
  };
  auto s2d = make_field("s2d",strip->get_2d_scalar_layout());
  auto v2d = make_field("v2d",strip->get_2d_vector_layout(vec_dim,"vec"));
  auto s3d = make_field("s3d",strip->get_3d_scalar_layout(LEV));
  auto v3d = make_field("v3d",strip->get_3d_vector_layout(LEV,vec_dim,"vec"));

  {
    const auto map = strip->get_unique_grid_gids().get_view<const gid_type*,Host>();
    auto h_s2d = s2d.get_view<Real*,Host>();
    auto h_v2d = v2d.get_view<Real**,Host>();
    auto h_s3d = s3d.get_view<Real**,Host>();
    auto h_v3d = v3d.get_view<Real***,Host>();
    for (int i=0; i<=n; ++i) {
      h_s2d(i) = value(map(i),0,0);
      for (int ic=0; ic<vec_dim; ++ic) {
        h_v2d(i,ic) = value(map(i),ic,0);
      }
      for (int il=0; il<num_levs; ++il) {
        h_s3d(i,il) = value(map(i),0,il);
        for (int ic=0; ic<vec_dim; ++ic) {
          h_v3d(i,ic,il) = value(map(i),ic,il);
        }
      }
    }
    for (auto f : {s2d,v2d,s3d,v3d}) {
      f.sync_to_dev();
    }
  }

  auto p_s2d = remap->register_field_from_src(s2d);
  auto p_v2d = remap->register_field_from_src(v2d);
  auto p_s3d = remap->register_field_from_src(s3d);
  auto p_v3d = remap->register_field_from_src(v3d);
  remap->registration_ends();

  for (auto f : {p_s2d,p_v2d,p_s3d,p_v3d}) {
    f.deep_copy(-1);
  }
  remap->remap_fwd();
  for (auto f : {p_s2d,p_v2d,p_s3d,p_v3d}) {
    f.sync_to_host();
  }

  const auto gids = unique->get_dofs_gids().get_view<const gid_type*,Host>();
  auto h_p_s2d = p_s2d.get_view<const Real*,Host>();
  auto h_p_v2d = p_v2d.get_view<const Real**,Host>();
  auto h_p_s3d = p_s3d.get_view<const Real**,Host>();
  auto h_p_v3d = p_v3d.get_view<const Real***,Host>();
  for (int icol=0; icol<n; ++icol) {
    CHECK (h_p_s2d(icol)==value(gids(icol),0,0));
    for (int ic=0; ic<vec_dim; ++ic) {
      CHECK (h_p_v2d(icol,ic)==value(gids(icol),ic,0));
    }
    for (int il=0; il<num_levs; ++il) {
      CHECK (h_p_s3d(icol,il)==value(gids(icol),0,il));
      for (int ic=0; ic<vec_dim; ++ic) {
        CHECK (h_p_v3d(icol,ic,il)==value(gids(icol),ic,il));
      }
    }
  }

  // Backward: set the unique fields, wipe the strip fields, and check that all the strip dofs
  // (including the one that is a replica of the first) get the value of their unique column
  for (auto f : {s2d,v2d,s3d,v3d}) {
    f.deep_copy(-1);
  }
  {
    auto h_s2d = p_s2d.get_view<Real*,Host>();
    auto h_v2d = p_v2d.get_view<Real**,Host>();
    auto h_s3d = p_s3d.get_view<Real**,Host>();
    auto h_v3d = p_v3d.get_view<Real***,Host>();
    for (int icol=0; icol<n; ++icol) {
      h_s2d(icol) = value(gids(icol),0,0);
      for (int ic=0; ic<vec_dim; ++ic) {
        h_v2d(icol,ic) = value(gids(icol),ic,0);
      }
      for (int il=0; il<num_levs; ++il) {
        h_s3d(icol,il) = value(gids(icol),0,il);
        for (int ic=0; ic<vec_dim; ++ic) {
          h_v3d(icol,ic,il) = value(gids(icol),ic,il);
        }
      }
    }
    for (auto f : {p_s2d,p_v2d,p_s3d,p_v3d}) {
      f.sync_to_dev();
    }
  }
  REQUIRE (remap->bwd_allowed());
  remap->remap_bwd();
  for (auto f : {s2d,v2d,s3d,v3d}) {
    f.sync_to_host();
  }
  {
    const auto map = strip->get_unique_grid_gids().get_view<const gid_type*,Host>();
    auto h_s2d = s2d.get_view<const Real*,Host>();
    auto h_v2d = v2d.get_view<const Real**,Host>();
    auto h_s3d = s3d.get_view<const Real**,Host>();
    auto h_v3d = v3d.get_view<const Real***,Host>();
    for (int i=0; i<=n; ++i) {
      CHECK (h_s2d(i)==value(map(i),0,0));
      for (int ic=0; ic<vec_dim; ++ic) {
        CHECK (h_v2d(i,ic)==value(map(i),ic,0));
      }
      for (int il=0; il<num_levs; ++il) {
        CHECK (h_s3d(i,il)==value(map(i),0,il));
        for (int ic=0; ic<vec_dim; ++ic) {
          CHECK (h_v3d(i,ic,il)==value(map(i),ic,il));
        }
      }
    }
  }
}

} // namespace scream
