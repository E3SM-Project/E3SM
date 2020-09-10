#include "share/grid/simple_grid.hpp"

namespace {
// Utility function, used in SimpleGrid ctor, to get around an issue with nvcc
void create_identity_map(scream::AbstractGrid::lid_to_idx_map_type view) {
  using RangePolicy = scream::AbstractGrid::kokkos_types::RangePolicy;
  Kokkos::parallel_for(RangePolicy(0,view.extent_int(0)),
                       KOKKOS_LAMBDA(const int i) {
    view(i,0) = i;
  });
}
} // anonymous namespace

namespace scream {

SimpleGrid::
SimpleGrid (const std::string& grid_name,
            const int num_global_cols,
            const int num_vertical_lev,
            const ekat::Comm& comm)
 : m_grid_name (grid_name)
 , m_num_vl    (num_vertical_lev)
{
  // Compute how many columns are owned by this rank
  const int num_procs = comm.size();

  m_num_my_cols = num_global_cols / num_procs;
  int remainder   = num_global_cols % num_procs;
  int dof_offset  = m_num_my_cols*comm.rank();
  if (comm.rank() < remainder) {
    ++m_num_my_cols;
    dof_offset += comm.rank();
  } else {
    dof_offset += remainder;
  }

  m_dofs_gids = decltype(m_dofs_gids)("phys dofs",m_num_my_cols);
  auto h_dofs_gids = Kokkos::create_mirror_view(m_dofs_gids);
  for (int i=0; i<m_num_my_cols; ++i) {
    h_dofs_gids(i) = dof_offset + i;
  }
  Kokkos::deep_copy(m_dofs_gids,h_dofs_gids);

  // The lid->idx map is the identity map.
  m_lid_to_idx = decltype(m_lid_to_idx)("lid to idx",m_num_my_cols,1);

  // Note: a KOKKOS_LAMBDA cannot be used in the c-tor, cause nvcc requires
  //       to be able to take the address of the enclosing function, and the
  //       c-tor of a class does not satisfy that. So wrap this in a utility
  //       free function, in an anonymous namespace 
  create_identity_map(m_lid_to_idx);
}

FieldLayout
SimpleGrid::get_native_dof_layout () const
{
  using namespace ShortFieldTagsNames;
  return FieldLayout({COL},{m_num_my_cols});
}

} // namespace scream
