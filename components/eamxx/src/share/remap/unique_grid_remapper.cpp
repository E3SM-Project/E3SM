#include "share/remap/unique_grid_remapper.hpp"

#include "share/util/eamxx_utils.hpp"

#include <ekat_pack_utils.hpp>
#include <ekat_assert.hpp>

namespace scream
{

KokkosTypes<DefaultDevice>::view_1d<int>
compute_unique_grid_dof_map (const AbstractGrid& grid,
                             const AbstractGrid& unique_grid)
{
  using gid_type = AbstractGrid::gid_type;

  const int num_unique_dofs = unique_grid.get_num_local_dofs();
  const int num_dofs        = grid.get_num_local_dofs();

  auto gids        = grid.get_unique_grid_gids().get_view<const gid_type*>();
  auto unique_gids = unique_grid.get_dofs_gids().get_view<const gid_type*>();

  auto policy = KokkosTypes<DefaultDevice>::RangePolicy(0,num_unique_dofs);
  KokkosTypes<DefaultDevice>::view_1d<int> u2g ("",num_unique_dofs);

  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idof){
    auto gid = unique_gids(idof);
    bool found = false;
    for (int i=0; i<num_dofs; ++i) {
      if (gids(i)==gid) {
        u2g(idof) = i;
        found = true;
        break;
      }
    }
    EKAT_KERNEL_ASSERT_MSG (found, "Error! Unique grid gid not found among the source grid dofs.\n");
    (void)found;
  });

  return u2g;
}

UniqueGridRemapper::
UniqueGridRemapper (const grid_ptr_type& grid)
 : AbstractRemapper(grid,grid->get_unique_grid())
{
  const auto unique_grid = m_tgt_grid;

  EKAT_REQUIRE_MSG (grid!=unique_grid,
      "Error! The source grid has no replicated dofs (it is its own unique grid). Use an IdentityRemapper.\n"
      " - grid name: " + grid->name() + "\n");
  EKAT_REQUIRE_MSG (grid->get_num_vertical_levels()==unique_grid->get_num_vertical_levels(),
      "Error! The source grid and its unique grid have a different number of vertical levels.\n"
      " - grid name       : " + grid->name() + "\n"
      " - unique grid name: " + unique_grid->name() + "\n");

  // Fields are addressed with K indices per dof. We only support K=1 or K=3 on the source
  // grid (e.g., (col) or (elem,gp,gp)), and K=1 on the target grid.
  auto get_num_idx = [](const AbstractGrid& g) {
    return g.get_lid_to_idx_map().get_header().get_identifier().get_layout().dim(1);
  };
  m_src_K = get_num_idx(*grid);
  EKAT_REQUIRE_MSG (m_src_K==1 or m_src_K==3,
      "Error! Only grids whose dofs are addressed with 1 or 3 indices are supported as source.\n"
      " - grid name: " + grid->name() + "\n"
      " - number of indices per dof: " + std::to_string(m_src_K) + "\n");
  EKAT_REQUIRE_MSG (get_num_idx(*unique_grid)==1,
      "Error! The target grid must have dofs addressed with a single index.\n"
      " - grid name: " + unique_grid->name() + "\n"
      " - number of indices per dof: " + std::to_string(get_num_idx(*unique_grid)) + "\n");

  m_bwd_allowed = false;

  m_num_unique_cols = unique_grid->get_num_local_dofs();
  m_src_lid2idx     = grid->get_lid_to_idx_map().get_view<const int**>();
  m_unique2src      = compute_unique_grid_dof_map(*grid,*unique_grid);
}

bool UniqueGridRemapper::
is_valid_src_layout (const FieldLayout& layout) const
{
  // We don't want fields with TimeLevel in it. Just subview the fields instead
  return AbstractRemapper::is_valid_src_layout(layout) and
         not ekat::contains(layout.tags(),FieldTag::TimeLevel);
}

void UniqueGridRemapper::
registration_ends_impl ()
{
  initialize_device_variables();
}

void UniqueGridRemapper::
initialize_device_variables ()
{
  m_src_repo.init(m_src_fields);
  m_tgt_repo.init(m_tgt_fields);

  m_layout              = decltype(m_layout)              ("layout", m_num_fields);
  m_pack_alloc_property = decltype(m_pack_alloc_property) ("pack_alloc_property", m_num_fields);
  m_num_levels          = decltype(m_num_levels)          ("num_physical_levels", m_num_fields);
  m_vec_dim             = decltype(m_vec_dim)             ("vec_dim", m_num_fields);

  auto h_layout              = Kokkos::create_mirror_view(m_layout);
  auto h_pack_alloc_property = Kokkos::create_mirror_view(m_pack_alloc_property);
  auto h_num_levels          = Kokkos::create_mirror_view(m_num_levels);
  auto h_vec_dim             = Kokkos::create_mirror_view(m_vec_dim);

  for (int i=0; i<m_num_fields; ++i) {
    const auto& sh = m_src_fields[i].get_header();
    const auto& th = m_tgt_fields[i].get_header();

    const auto& sl = sh.get_identifier().get_layout();

    auto lt = sl.type();
    h_layout(i) = etoi(lt);

    const bool is_field_3d = lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D;
    h_num_levels(i) = is_field_3d ? sl.dims().back() : -1;
    h_vec_dim(i)    = sl.is_vector_layout() ? sl.get_vector_dim() : -1;

    const auto& sap = sh.get_alloc_properties();
    const auto& tap = th.get_alloc_properties();
    if (is_field_3d &&
        sap.template is_compatible<pack_type>() &&
        tap.template is_compatible<pack_type>()) {
      h_pack_alloc_property(i) = AllocPropType::PackAlloc;
    } else {
      h_pack_alloc_property(i) = AllocPropType::RealAlloc;
    }
  }
  Kokkos::deep_copy(m_layout,              h_layout             );
  Kokkos::deep_copy(m_pack_alloc_property, h_pack_alloc_property);
  Kokkos::deep_copy(m_num_levels,          h_num_levels         );
  Kokkos::deep_copy(m_vec_dim,             h_vec_dim            );
}

void UniqueGridRemapper::
remap_fwd_impl ()
{
  // Check if we need to update the views for subfields
  m_src_repo.update_subfields(m_src_fields);
  m_tgt_repo.update_subfields(m_tgt_fields);

  // TeamPolicy over m_num_unique_cols*m_num_fields
  using TeamPolicy = typename KT::TeamTagPolicy<RemapTag>;
  const TeamPolicy policy(m_num_fields*m_num_unique_cols,Kokkos::AUTO);
  Kokkos::parallel_for(policy, *this);
  Kokkos::fence();
}

// NOTE: in the kernels below, the source field is addressed with the dof 'lid' (if m_src_K==1),
//       or with the three indices (e,g1,g2) in m_src_lid2idx (if m_src_K==3). In the latter case,
//       the indices are in the same order as the horizontal tags in the grid's 2d scalar layout.
//       Any other dim in the layout is located by the layout type. E.g., for SE grids,
//       a vector dim sits between the first and second index.

template <typename MT>
KOKKOS_FUNCTION
void UniqueGridRemapper::
local_remap_2d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int icol = rank / m_num_fields;
  const int lid  = m_unique2src(icol);
  const bool K3  = m_src_K==3;

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto u = m_tgt_repo.views[i].v1d;
      if (K3) {
        auto s = m_src_repo.cviews[i].v3d;
        u(icol) = s(m_src_lid2idx(lid,0),m_src_lid2idx(lid,1),m_src_lid2idx(lid,2));
      } else {
        auto s = m_src_repo.cviews[i].v1d;
        u(icol) = s(lid);
      }
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      auto u = m_tgt_repo.views[i].v2d;
      const int vec_dim = m_vec_dim(i);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim);
      if (K3) {
        auto s = m_src_repo.cviews[i].v4d;
        const int e = m_src_lid2idx(lid,0), g1 = m_src_lid2idx(lid,1), g2 = m_src_lid2idx(lid,2);
        Kokkos::parallel_for(tr, [&] (const int idim) {
          u(icol,idim) = s(e,idim,g1,g2);
        });
      } else {
        auto s = m_src_repo.cviews[i].v2d;
        Kokkos::parallel_for(tr, [&] (const int idim) {
          u(icol,idim) = s(lid,idim);
        });
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template <typename ScalarT, typename MT>
KOKKOS_FUNCTION
void UniqueGridRemapper::
local_remap_3d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int icol = rank / m_num_fields;
  const int lid  = m_unique2src(icol);
  const bool K3  = m_src_K==3;

  constexpr int PackSize = sizeof(ScalarT) / sizeof(Real);
  using PI = ekat::PackInfo<PackSize>;
  const int num_packs = PI::num_packs(m_num_levels(i));

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar3D):
    {
      auto u = pack_view<      ScalarT>(m_tgt_repo.views[i].v2d);

      const auto tr = Kokkos::TeamVectorRange(team, num_packs);
      if (K3) {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v4d);
        const int e = m_src_lid2idx(lid,0), g1 = m_src_lid2idx(lid,1), g2 = m_src_lid2idx(lid,2);
        Kokkos::parallel_for(tr, [&] (const int ilev) {
          u(icol,ilev) = s(e,g1,g2,ilev);
        });
      } else {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v2d);
        Kokkos::parallel_for(tr, [&] (const int ilev) {
          u(icol,ilev) = s(lid,ilev);
        });
      }
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto u = pack_view<      ScalarT>(m_tgt_repo.views[i].v3d);
      const int vec_dim = m_vec_dim(i);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim*num_packs);
      if (K3) {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v5d);
        const int e = m_src_lid2idx(lid,0), g1 = m_src_lid2idx(lid,1), g2 = m_src_lid2idx(lid,2);
        Kokkos::parallel_for(tr, [&] (const int idx) {
          const int idim = idx / num_packs;
          const int ilev = idx % num_packs;
          u(icol,idim,ilev) = s(e,idim,g1,g2,ilev);
        });
      } else {
        auto s = pack_view<const ScalarT>(m_src_repo.cviews[i].v3d);
        Kokkos::parallel_for(tr, [&] (const int idx) {
          const int idim = idx / num_packs;
          const int ilev = idx % num_packs;
          u(icol,idim,ilev) = s(lid,idim,ilev);
        });
      }
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename MT>
KOKKOS_INLINE_FUNCTION
void UniqueGridRemapper::
operator()(const RemapTag&, const MT& team) const
{
  const int rank = team.league_rank();
  const int i = rank % m_num_fields;

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    case etoi(LayoutType::Vector2D):
      local_remap_2d(team);
      break;
    case etoi(LayoutType::Scalar3D):
    case etoi(LayoutType::Vector3D):
      if (m_pack_alloc_property(i) == AllocPropType::PackAlloc) {
        local_remap_3d<pack_type>(team);
      } else {
        local_remap_3d<Real>(team);
      }
      break;
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

} // namespace scream
