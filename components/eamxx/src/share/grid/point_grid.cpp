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

PointGrid::
PointGrid (const std::string& grid_name,
           const int          num_my_cols,
           const int          num_global_cols,
           const int          num_vertical_levels,
           const ekat::Comm&  comm)
 : AbstractGrid(grid_name,GridType::Point,num_my_cols,num_global_cols,num_vertical_levels,comm)
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

  return FieldLayout({COL},{get_num_local_dofs()}).rename_dims(m_special_tag_names);
}

FieldLayout
PointGrid::get_2d_vector_layout (const int vector_dim, const std::string& vec_dim_name) const
{
  using namespace ShortFieldTagsNames;

  FieldLayout fl({COL,CMP},{get_num_local_dofs(),vector_dim});
  fl.rename_dim(1,vec_dim_name);
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout
PointGrid::get_2d_tensor_layout (const std::vector<int>& cmp_dims,
                                 const std::vector<std::string>& cmp_names) const
{
  EKAT_REQUIRE_MSG (cmp_names.size()==cmp_dims.size(),
      "[PointGrid::get_2d_tensor_layout] Input vector dimensions mismatch.\n"
      "  - grid name: " + name() + "\n"
      "  - cmp_names: " + ekat::join(cmp_names,",") + "\n"
      "  - cmp_dims : " + ekat::join(cmp_dims,",") + "\n");
  using namespace ShortFieldTagsNames;

  FieldLayout fl;

  fl.append_dim(COL, get_num_local_dofs());
  for (size_t i=0; i<cmp_dims.size(); ++i) {
    fl.append_dim(CMP,cmp_dims[i],cmp_names[i]);
  }

  return fl.rename_dims(m_special_tag_names);
}

FieldLayout
PointGrid::get_3d_scalar_layout (const bool midpoints) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  return FieldLayout({COL,VL},{get_num_local_dofs(),nvl}).rename_dims(m_special_tag_names);
}

FieldLayout
PointGrid::get_3d_vector_layout (const bool midpoints, const int vector_dim,
                                 const std::string& vec_dim_name) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  FieldLayout fl({COL,CMP,VL},{get_num_local_dofs(),vector_dim,nvl});
  fl.rename_dim(1,vec_dim_name);
  return fl.rename_dims(m_special_tag_names);
}

FieldLayout
PointGrid::get_3d_tensor_layout (const bool midpoints,
                                 const std::vector<int>& cmp_dims,
                                 const std::vector<std::string>& cmp_names) const
{
  EKAT_REQUIRE_MSG (cmp_names.size()==cmp_dims.size(),
      "[PointGrid::get_2d_tensor_layout] Input vector dimensions mismatch.\n"
      "  - grid name: " + name() + "\n"
      "  - cmp_names: " + ekat::join(cmp_names,",") + "\n"
      "  - cmp_dims : " + ekat::join(cmp_dims,",") + "\n");
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);
  auto VL = midpoints ? LEV : ILEV;

  FieldLayout fl;

  fl.append_dim(COL, get_num_local_dofs());
  for (size_t i=0; i<cmp_dims.size(); ++i) {
    fl.append_dim(CMP,cmp_dims[i],cmp_names[i]);
  }
  fl.append_dim(VL,nvl);

  return fl.rename_dims(m_special_tag_names);
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
