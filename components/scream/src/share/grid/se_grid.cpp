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
 , m_num_local_dofs (0)
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
  std::vector<FieldTag> tags;
  std::vector<int> dims;
  switch (m_type) {
    case GridType::SE_CellBased:
      tags = {EL,GP,GP};
      dims = {m_num_local_elem,m_num_gp,m_num_gp};
      break;
    case GridType::SE_NodeBased:
      tags = {COL};
      dims = {m_num_local_dofs};
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unexpected grid type. Please, contact developers.\n");
  }

  return FieldLayout(tags,dims);
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

  const int expected_dim_1 = m_type==GridType::SE_CellBased ? 3 : 1;

  EKAT_REQUIRE_MSG (
    lid_to_idx.extent_int(1)==expected_dim_1,
    "Error! Invalid extent(1) for lid_to_idx input:\n"
      "   expected: " + std::to_string(expected_dim_1) + "\n"
      "   provided: " + std::to_string(lid_to_idx.extent_int(1)) + "\n");

  m_num_local_dofs = dofs.extent_int(0);
  if (m_type==GridType::SE_CellBased) {
    EKAT_REQUIRE_MSG (
      m_num_local_dofs==(m_num_local_elem*m_num_gp*m_num_gp),
      "Error! Input dofs view has the wrong size:\n"
      "   expected: " + std::to_string(m_num_local_elem*m_num_gp*m_num_gp) + "\n"
      "   provided: " + std::to_string(dofs.extent_int(0)) + "\n");
  }
  m_dofs_gids  = dofs;
  m_lid_to_idx = lid_to_idx;
}

void SEGrid::create_elgpgp_to_lid_map ()
{
  // Create idx->gid map, and init with invalid entries (-1)
  // For CellBased, all -1 should be overwritten, but NodeBased, in parallel,
  // might have some -1 left (useful to detect access to off-rank entries)
  m_elgpgp_to_lid = decltype(m_elgpgp_to_lid)("elgpgp to lid",m_num_local_elem,m_num_gp,m_num_gp);

  auto lid_to_idx = m_lid_to_idx;
  Kokkos::deep_copy(m_elgpgp_to_lid,-1);

  // Loop over lid_to_elgpgp, and fill available entries from lid_to_elgpgp
  auto elgpgp_to_lid = m_elgpgp_to_lid;
  Kokkos::parallel_for(kokkos_types::RangePolicy(0,m_dofs_gids.extent_int(0)),
                       KOKKOS_LAMBDA(const int i) {
    const int ie  = lid_to_idx(i,0);
    const int igp = lid_to_idx(i,1);
    const int jgp = lid_to_idx(i,2);

    elgpgp_to_lid(ie,igp,jgp) = i;
  });
}

} // namespace scream
