#include "share/grid/se_grid.hpp"

#include "ekat/kokkos//ekat_subview_utils.hpp"

namespace scream {

SEGrid::
SEGrid (const std::string& grid_name,
        const int num_my_elements,
        const int num_gauss_pts,
        const int num_vertical_levels,
        const ekat::Comm& comm)
 : AbstractGrid (grid_name,GridType::SE,num_my_elements*num_gauss_pts*num_gauss_pts,num_vertical_levels,comm)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (num_my_elements>=0, "Error! Number of local elements must be non-negative.\n");
  EKAT_REQUIRE_MSG (num_gauss_pts>=2, "Error! Number of gauss points must be at least 2.\n");
  EKAT_REQUIRE_MSG (num_vertical_levels>=2, "Error! Number of vertical levels must be at least 2.\n");

  m_num_local_elem = num_my_elements;
  m_num_gp         = num_gauss_pts;
  get_comm().all_reduce(&m_num_local_elem,&m_num_global_elem,1,MPI_SUM);
}

FieldLayout
SEGrid::get_2d_scalar_layout () const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({EL,GP,GP},{m_num_local_elem,m_num_gp,m_num_gp});
}

FieldLayout
SEGrid::get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({EL,vector_tag,GP,GP},{m_num_local_elem,vector_dim,m_num_gp,m_num_gp});
}

FieldLayout
SEGrid::get_3d_scalar_layout (const bool midpoints) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({EL,GP,GP,VL},{m_num_local_elem,m_num_gp,m_num_gp,nvl});
}

FieldLayout
SEGrid::get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({EL,vector_tag,GP,GP,VL},{m_num_local_elem,vector_dim,m_num_gp,m_num_gp,nvl});
}

auto SEGrid::get_cg_dofs_gids () const -> const dofs_list_type&
{
  EKAT_REQUIRE_MSG (m_cg_dofs_set, "Error! CG dofs not yet set.\n");
  return m_cg_dofs_gids;
}

void SEGrid::set_cg_dofs (const dofs_list_type& cg_dofs)
{
  EKAT_REQUIRE_MSG (not m_cg_dofs_set, "Error! CG dofs were already set.\n");
  m_cg_dofs_gids = cg_dofs;
  m_cg_dofs_set = true;
}

std::shared_ptr<AbstractGrid> SEGrid::clone (const std::string& clone_name, const bool shallow) const
{
  auto grid = std::make_shared<SEGrid>(clone_name,m_num_local_elem,m_num_gp,get_num_vertical_levels(),get_comm());

  grid->copy_views(*this,shallow);

  if (m_cg_dofs_set) {
    if (shallow) {
      grid->set_cg_dofs(m_cg_dofs_gids);
    } else {
      decltype(m_cg_dofs_gids) cg_dofs ("",m_cg_dofs_gids.size());
      Kokkos::deep_copy(cg_dofs,m_cg_dofs_gids);
      grid->set_cg_dofs(cg_dofs);
    }
  }

  return grid;
}

bool SEGrid::valid_dofs_list (const dofs_list_type& /*dofs_gids*/) const
{
  return is_unique();
}

bool SEGrid::valid_lid_to_idx_map (const lid_to_idx_map_type& /*lid_to_idx*/) const
{
  auto l2i = get_lid_to_idx_map();
  auto h_l2i = Kokkos::create_mirror_view(l2i);
  Kokkos::deep_copy(h_l2i,l2i);
  for (int idof=0; idof<h_l2i.extent_int(0); ++idof) {
    auto elgpgp = ekat::subview(h_l2i,idof);
    const int el = elgpgp(0);
    const int ip = elgpgp(1);
    const int jp = elgpgp(2);

    if (el<0 || el>=m_num_local_elem) return false;
    if (ip<0 || ip>=m_num_gp) return false;
    if (jp<0 || jp>=m_num_gp) return false;
  }
  return true;
}

} // namespace scream
