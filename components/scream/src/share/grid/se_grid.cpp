#include "share/grid/se_grid.hpp"

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
}

SEGrid::
SEGrid (const std::string& grid_name,
        const int num_my_elements,
        const int num_gauss_pts,
        const int num_vertical_levels,
        const std::shared_ptr<const AbstractGrid>& unique_grid,
        const ekat::Comm& comm)
 : SEGrid(grid_name,num_my_elements,num_gauss_pts,num_vertical_levels,comm)
{
  set_unique_grid (unique_grid);
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

} // namespace scream
