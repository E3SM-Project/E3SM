#include "share/grid/point_grid.hpp"

#include <numeric>

namespace scream {

PointGrid::
PointGrid (const std::string&     grid_name,
           const dofs_list_type&  dofs_gids,
           const int              num_vertical_lev)
 : m_grid_name   (grid_name)
 , m_num_my_cols (dofs_gids.size())
 , m_num_vl      (num_vertical_lev)
 , m_dofs_gids   (dofs_gids)
{
  // The lid->idx map is the identity map.
  m_lid_to_idx = decltype(m_lid_to_idx)("lid to idx",m_num_my_cols,1);

  auto h_lid_to_idx = Kokkos::create_mirror_view(m_lid_to_idx);
  std::iota(h_lid_to_idx.data(),h_lid_to_idx.data()+m_num_my_cols,0);
  Kokkos::deep_copy(m_lid_to_idx,h_lid_to_idx);
}

FieldLayout
PointGrid::get_native_dof_layout () const
{
  using namespace ShortFieldTagsNames;
  return FieldLayout({COL},{m_num_my_cols});
}

PointGrid
create_point_grid (const std::string& grid_name,
                   const int num_global_cols,
                   const int num_vertical_lev,
                   const ekat::Comm& comm)
{
  // Compute how many columns are owned by this rank
  const int num_procs = comm.size();

  auto num_my_cols = num_global_cols / num_procs;
  int remainder   = num_global_cols % num_procs;
  int dof_offset  = num_my_cols*comm.rank();
  if (comm.rank() < remainder) {
    ++num_my_cols;
    dof_offset += comm.rank();
  } else {
    dof_offset += remainder;
  }

  PointGrid::dofs_list_type dofs_gids ("phys dofs",num_my_cols);
  auto h_dofs_gids = Kokkos::create_mirror_view(dofs_gids);
  for (int i=0; i<num_my_cols; ++i) {
    h_dofs_gids(i) = dof_offset + i;
  }
  Kokkos::deep_copy(dofs_gids,h_dofs_gids);

  return PointGrid(grid_name,dofs_gids,num_vertical_lev);
}

} // namespace scream
