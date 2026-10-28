#include "share/remap/dg_to_cg_remapper.hpp"

#include "share/grid/se_grid.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat_pack_utils.hpp>
#include <ekat_assert.hpp>

namespace scream
{

KokkosTypes<DefaultDevice>::view_1d<int>
compute_cg_to_se_dof_map (const AbstractGrid& se_grid,
                          const AbstractGrid& cg_grid)
{
  using gid_type = AbstractGrid::gid_type;

  const int num_cg_dofs = cg_grid.get_num_local_dofs();
  const int num_se_dofs = se_grid.get_num_local_dofs();

  auto se_gids = se_grid.get_lid_to_unique_gid_map().get_view<const gid_type*>();
  auto cg_gids = cg_grid.get_dofs_gids().get_view<const gid_type*>();

  auto policy = KokkosTypes<DefaultDevice>::RangePolicy(0,num_cg_dofs);
  KokkosTypes<DefaultDevice>::view_1d<int> p2d ("",num_cg_dofs);

  Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const int idof){
    auto gid = cg_gids(idof);
    bool found = false;
    for (int i=0; i<num_se_dofs; ++i) {
      if (se_gids(i)==gid) {
        p2d(idof) = i;
        found = true;
        break;
      }
    }
    EKAT_KERNEL_ASSERT_MSG (found, "Error! CG grid gid not found in the SE grid.\n");
    (void)found;
  });

  return p2d;
}

DGToCGRemapper::
DGToCGRemapper (const grid_ptr_type& se_grid,
                const grid_ptr_type& cg_grid)
 : AbstractRemapper(se_grid,cg_grid)
{
  EKAT_REQUIRE_MSG(se_grid->type()==GridType::SE,    "Error! Input source grid is not a SE grid.\n");
  EKAT_REQUIRE_MSG(cg_grid->type()==GridType::Point, "Error! Input target grid is not a Point grid.\n");

  m_bwd_allowed = false;

  m_num_cg_cols = cg_grid->get_num_local_dofs();
  m_lid2elgp    = se_grid->get_lid_to_idx_map().get_view<const int**>();
  m_cg2se         = compute_cg_to_se_dof_map(*se_grid,*cg_grid);
}

bool DGToCGRemapper::
is_valid_src_layout (const FieldLayout& layout) const
{
  // We don't want fields with TimeLevel in it. Just subview the fields instead
  return AbstractRemapper::is_valid_src_layout(layout) and
         not ekat::contains(layout.tags(),FieldTag::TimeLevel);
}

void DGToCGRemapper::
registration_ends_impl ()
{
  initialize_device_variables();
}

void DGToCGRemapper::
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

void DGToCGRemapper::
remap_fwd_impl ()
{
  // Check if we need to update the views for subfields
  m_src_repo.update_subfields(m_src_fields);
  m_tgt_repo.update_subfields(m_tgt_fields);

  // TeamPolicy over m_num_cg_cols*m_num_fields
  using TeamPolicy = typename KT::TeamTagPolicy<RemapTag>;
  const TeamPolicy policy(m_num_fields*m_num_cg_cols,Kokkos::AUTO);
  Kokkos::parallel_for(policy, *this);
  Kokkos::fence();
}

template <typename MT>
KOKKOS_FUNCTION
void DGToCGRemapper::
local_remap_2d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int icol = rank / m_num_fields;
  const auto& elgp = Kokkos::subview(m_lid2elgp,m_cg2se(icol),Kokkos::ALL());

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar2D):
    {
      auto cg = m_tgt_repo.views[i].v1d;
      auto se = m_src_repo.cviews[i].v3d;

      cg(icol) = se(elgp[0],elgp[1],elgp[2]);
      break;
    }
    case etoi(LayoutType::Vector2D):
    {
      auto cg = m_tgt_repo.views[i].v2d;
      auto se = m_src_repo.cviews[i].v4d;
      const int vec_dim = m_vec_dim(i);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim);
      const auto f = [&] (const int idim) {
        cg(icol,idim) = se(elgp[0],idim,elgp[1],elgp[2]);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template <typename ScalarT, typename MT>
KOKKOS_FUNCTION
void DGToCGRemapper::
local_remap_3d (const MT& team) const
{
  const int rank = team.league_rank();
  const int i    = rank % m_num_fields;
  const int icol = rank / m_num_fields;
  const auto& elgp = Kokkos::subview(m_lid2elgp,m_cg2se(icol),Kokkos::ALL());

  constexpr int PackSize = sizeof(ScalarT) / sizeof(Real);
  using PI = ekat::PackInfo<PackSize>;
  const int num_packs = PI::num_packs(m_num_levels(i));

  switch (m_layout(i)) {
    case etoi(LayoutType::Scalar3D):
    {
      auto cg = pack_view<      ScalarT>(m_tgt_repo.views[i].v2d);
      auto se = pack_view<const ScalarT>(m_src_repo.cviews[i].v4d);

      const auto tr = Kokkos::TeamVectorRange(team, num_packs);
      const auto f = [&] (const int ilev) {
        cg(icol,ilev) = se(elgp[0],elgp[1],elgp[2],ilev);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    case etoi(LayoutType::Vector3D):
    {
      auto cg = pack_view<      ScalarT>(m_tgt_repo.views[i].v3d);
      auto se = pack_view<const ScalarT>(m_src_repo.cviews[i].v5d);
      const int vec_dim = m_vec_dim(i);

      const auto tr = Kokkos::TeamVectorRange(team, vec_dim*num_packs);
      const auto f = [&] (const int idx) {
        const int idim = idx / num_packs;
        const int ilev = idx % num_packs;
        cg(icol,idim,ilev) = se(elgp[0],idim,elgp[1],elgp[2],ilev);
      };
      Kokkos::parallel_for(tr, f);
      break;
    }
    default:
      EKAT_KERNEL_ERROR_MSG("Error! Unhandled case in switch statement.\n");
  }
}

template<typename MT>
KOKKOS_INLINE_FUNCTION
void DGToCGRemapper::
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
