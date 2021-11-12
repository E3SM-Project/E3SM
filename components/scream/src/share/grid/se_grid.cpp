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
}

void SEGrid::
check_geo_data (const std::string& name, const geo_view_type& /* data */) const {
  // Sanity checks
  EKAT_REQUIRE_MSG (name=="lat" || name=="lon" || name=="area",
                    "Error! Point grid does not support geometry data '" + name + "'.\n");

  // TODO: check actual values?
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

void SEGrid::check_lid_to_idx_map () const
{
  auto l2i = get_lid_to_idx_map();
  auto h_l2i = Kokkos::create_mirror_view(l2i);
  Kokkos::deep_copy(h_l2i,l2i);
  for (int idof=0; idof<h_l2i.extent_int(0); ++idof) {
    auto elgpgp = ekat::subview(h_l2i,idof);
    const int el = elgpgp(0);
    const int ip = elgpgp(1);
    const int jp = elgpgp(2);

    EKAT_REQUIRE_MSG (el>=0 && el<m_num_local_elem,
        "Error! Element index out of bounds.\n");
    EKAT_REQUIRE_MSG (ip>=0 && ip<m_num_gp,
        "Error! Gauss point index i out of bounds.\n");
    EKAT_REQUIRE_MSG (jp>=0 && jp<m_num_gp,
        "Error! Gauss point index j out of bounds.\n");
  }
}

} // namespace scream
