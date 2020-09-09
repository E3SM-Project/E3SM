#include <share/grid/se_grid.hpp>

namespace scream {

SEGrid::
SEGrid (const std::string& grid_name,
        const GridType type,
        const int num_local_elements,
        const int num_gp,
        const int num_vl)
 : m_grid_name      (grid_name)
 , m_type           (type)
 , m_num_local_elem (num_local_elements)
 , m_num_gp         (num_gp)
 , m_num_vl         (num_vl)
 , m_num_local_dofs (-1)
{
  // Sanity check
  EKAT_REQUIRE_MSG (
    m_type==GridType::SE_CellBased || m_type==GridType::SE_NodeBased,
    "Error! Invalid grid type: " + e2str(m_type) + "\n");
}

FieldLayout
SEGrid::get_native_dof_layout () const
{
  using namespace ShortFieldTagsNames;
  switch (m_type) {
    case GridType::SE_CellBased:
      return FieldLayout({EL,GP,GP},{m_num_local_elem,m_num_gp,m_num_gp});
    case GridType::SE_NodeBased:
      return FieldLayout({COL},{m_num_local_dofs});
    default:
      EKAT_ERROR_MSG ("Error! Unexpected grid type. Please, contact developers.\n");
  }

  return FieldLayout({});
}

void SEGrid::
set_dofs (const dofs_list_type&      dofs,
          const lid_to_idx_map_type& lid_to_elgpgp)
{
  EKAT_REQUIRE_MSG (
    dofs.extent_int(0)==lid_to_elgpgp.extent_int(0),
    "Error! Input views dimension do not match:\n"
    "         num gids: " + std::to_string(dofs.extent_int(0)) + "\n"
    "         num lids: " + std::to_string(lid_to_elgpgp.extent_int(0)) + "\n");

  m_num_local_dofs = dofs.extent_int(0);
  if (m_type==GridType::SE_CellBased) {
    EKAT_REQUIRE_MSG (
      m_num_local_dofs==(m_num_local_elem*m_num_gp*m_num_gp),
      "Error! Input dofs view has the wrong size:\n"
      "   expected: " + std::to_string(m_num_local_elem*m_num_gp*m_num_gp) + "\n"
      "   provided: " + std::to_string(dofs.extent_int(0)) + "\n");
  }
  m_dofs_gids     = dofs;
  m_lid_to_elgpgp = lid_to_elgpgp;

  // Create idx->gid map, and init with invalid entries (-1)
  m_elgpgp_to_gid = decltype(m_elgpgp_to_gid)("elgpgp to gid",m_num_local_elem,m_num_gp,m_num_gp);
  Kokkos::deep_copy(m_elgpgp_to_gid,-1);

  // Loop over lid_to_elgpgp, and fill available entries from lid_to_elgpgp
  auto elgpgp_to_gid = m_elgpgp_to_gid;
  Kokkos::parallel_for(kokkos_types::RangePolicy(0,dofs.extent_int(0)),
                       KOKKOS_LAMBDA(const int i) {
    const int ie  = lid_to_elgpgp(i,0);
    const int igp = lid_to_elgpgp(i,1);
    const int jgp = lid_to_elgpgp(i,2);

    elgpgp_to_gid(ie,igp,jgp) = i;
  });
}

} // namespace scream
