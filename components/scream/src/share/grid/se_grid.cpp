#include "share/grid/se_grid.hpp"

namespace scream {

SEGrid::
SEGrid (const std::string& grid_name,
        const int num_global_elements,
        const int num_my_elements,
        const int num_gauss_pts,
        const int num_vertical_levels)
 : AbstractGrid (GridType::SE, grid_name)
 , m_num_local_elem (num_my_elements)
 , m_num_gp         (num_gauss_pts)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (
    num_global_elements>=num_my_elements,
    "Error! The number of local element is larger than the global one.\n"
    "       Did you call the constructor with the two swapped?\n");

  m_num_local_dofs  = m_num_local_elem*m_num_gp*m_num_gp;
  m_num_global_dofs = num_global_elements*m_num_gp*m_num_gp;
  m_num_vert_levs   = num_vertical_levels;
}

void SEGrid::
set_dofs (const dofs_list_type&      dofs,
          const lid_to_idx_map_type& lid_to_idx)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (
    dofs.extent_int(0)==lid_to_idx.extent_int(0),
    "Error! Input views dimension do not match:\n"
    "         num gids: " + std::to_string(dofs.extent_int(0)) + "\n"
    "         num lids: " + std::to_string(lid_to_idx.extent_int(0)) + "\n");

  EKAT_REQUIRE_MSG (
    dofs.extent_int(0)==m_num_local_dofs,
    "Error! Input views have the wrong size:\n"
    "   expected: " + std::to_string(m_num_local_dofs) + "\n"
    "   provided: " + std::to_string(dofs.extent_int(0)) + "\n");

  EKAT_REQUIRE_MSG (
    lid_to_idx.extent_int(1)==3,
    "Error! Invalid extent(1) for lid_to_idx input:\n"
      "   expected: 3\n"
      "   provided: " + std::to_string(lid_to_idx.extent_int(1)) + "\n");

  m_dofs_gids  = dofs;
  m_lid_to_idx = lid_to_idx;
}

void SEGrid::
set_geometry_data (const std::string& name, const geo_view_type& data) {
  // Sanity checks
  EKAT_REQUIRE_MSG (data.extent_int(0)==m_num_local_dofs,
                    "Error! Input geometry data has wrong dimensions.\n");
  EKAT_REQUIRE_MSG (name=="lat" || name=="lon" || name=="area",
                    "Error! Point grid does not support geometry data '" + name + "'.\n");

  m_geo_views[name] = data;
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

  return FieldLayout({EL,GP,GP,VL},{m_num_local_elem,m_num_gp,m_num_gp,nvl});
}

FieldLayout
SEGrid::get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const
{
  using namespace ShortFieldTagsNames;

  int nvl = this->get_num_vertical_levels() + (midpoints ? 0 : 1);

  return FieldLayout({EL,vector_tag,GP,GP,VL},{m_num_local_elem,vector_dim,m_num_gp,m_num_gp,nvl});
}

} // namespace scream
