#include <share/grid/se_grid.hpp>

namespace scream {

SEGrid::
SEGrid (const std::string& grid_name,
        const int num_local_elements,
        const int num_gauss_pts,
        const int num_vertical_levels)
 : m_grid_name      (grid_name)
 , m_num_local_elem (num_local_elements)
 , m_num_gp         (num_gauss_pts)
 , m_num_vl         (num_vertical_levels)
{
   m_num_local_dofs = m_num_local_elem*m_num_gp*m_num_gp;
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

} // namespace scream
