#include "share/grid/point_grid.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/share/physics_constants.hpp"

#include <numeric>

namespace scream {

PointGrid::
PointGrid (const std::string& grid_name,
           const int          num_my_cols,
           const int          num_vertical_levels,
           const ekat::Comm&  comm)
 : AbstractGrid(grid_name,GridType::Point,num_my_cols,num_vertical_levels,comm)
{
  create_dof_fields (get_2d_scalar_layout().rank());

  // The lid->idx map is the identity map.
  auto lid2idx = get_lid_to_idx_map();
  auto h_lid_to_idx = lid2idx.get_view<int**,Host>();
  std::iota(h_lid_to_idx.data(),h_lid_to_idx.data()+get_num_local_dofs(),0);
  lid2idx.sync_to_dev();
}

FieldLayout
PointGrid::get_2d_scalar_layout () const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({COL},{get_num_local_dofs()});
}

FieldLayout
PointGrid::get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({COL,vector_tag},{get_num_local_dofs(),vector_dim});
}

FieldLayout
PointGrid::get_3d_scalar_layout (const bool midpoints) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({COL,VL},{get_num_local_dofs(),nvl});
}

FieldLayout
PointGrid::get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({COL,vector_tag,VL},{get_num_local_dofs(),vector_dim,nvl});
}

std::shared_ptr<AbstractGrid>
PointGrid::clone (const std::string& clone_name,
                  const bool shallow) const
{
  auto grid = std::make_shared<PointGrid> (clone_name,get_num_local_dofs(),get_num_vertical_levels(),get_comm());
  grid->copy_data(*this,shallow);
  return grid;
}

std::shared_ptr<PointGrid>
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

  auto grid = std::make_shared<PointGrid>(grid_name,num_my_cols,num_vertical_lev,comm);
  grid->setSelfPointer(grid);

  auto dofs_gids = grid->get_dofs_gids();
  auto h_dofs_gids = dofs_gids.get_view<AbstractGrid::gid_type*,Host>();
  std::iota(h_dofs_gids.data(),h_dofs_gids.data()+num_my_cols,dof_offset);
  dofs_gids.sync_to_dev();

  return grid;
}

} // namespace scream
