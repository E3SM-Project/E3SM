#ifndef SCREAM_SE_GRID_HPP
#define SCREAM_SE_GRID_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"
#include "ekat/ekat_assert.hpp"

namespace scream
{

class SEGrid : public AbstractGrid
{
public:

  using elgpgp_to_lid_map_type = kokkos_types::view<int***>; // (elem, igp, jgp) -> int

  SEGrid (const std::string& grid_name,
          const GridType type,
          const int num_local_elements,
          const int num_gauss_pts,
          const int num_vertical_levels);

  virtual ~SEGrid () = default;

  // Grid description utilities
  GridType type () const override { return m_type; }
  const std::string& name () const override { return m_grid_name; }

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a 2d structured grid, this could be a set of 2 indices.
  FieldLayout get_native_dof_layout () const override;

  int get_num_vertical_levels () const override { return m_num_vl; }

  // Dofs gids utilities
  int get_num_local_dofs () const override { return m_num_local_dofs; }
  const dofs_list_type& get_dofs_gids () const override { return m_dofs_gids; }
  lid_to_idx_map_type get_lid_to_idx_map () const override { return m_lid_to_idx; }

  // Methods specific to SEGrid
  void create_elgpgp_to_lid_map ();
  elgpgp_to_lid_map_type get_idx_to_lid_map () const { return m_elgpgp_to_lid; }
  void set_dofs (const dofs_list_type&      dofs,
                 const lid_to_idx_map_type& lid_to_elgpgp);

protected:

  // Description variables
  const std::string         m_grid_name;
  const GridType            m_type;

  // Counters
  int                       m_num_local_elem;
  int                       m_num_gp;
  int                       m_num_vl;
  int                       m_num_local_dofs;

  // Dofs info variables
  dofs_list_type            m_dofs_gids;
  lid_to_idx_map_type       m_lid_to_idx;
  elgpgp_to_lid_map_type    m_elgpgp_to_lid;
};

} // namespace scream

#endif // SCREAM_SE_GRID_HPP
