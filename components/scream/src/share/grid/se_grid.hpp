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

  SEGrid (const std::string& grid_name,
          const int num_local_elements,
          const int num_gauss_pts,
          const int num_vertical_levels);

  virtual ~SEGrid () = default;

  // Grid description utilities
  GridType type () const override { return GridType::SE; }
  const std::string& name () const override { return m_grid_name; }

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a 2d structured grid, this could be a set of 2 indices.
  FieldLayout get_native_dof_layout () const override;

  int get_num_vertical_levels () const override { return m_num_vl; }

  // Dofs gids utilities
  int get_num_local_dofs () const override { return m_num_local_dofs; }
  const dofs_list_type& get_dofs_gids () const override { return m_dofs_gids; }
  lid_to_idx_map_type get_lid_to_idx_map () const override { return m_lid_to_idx; }

  // Set the dofs gids as well as their coordinates (ie,igp,jgp) in the grid
  void set_dofs (const dofs_list_type&      dofs,
                 const lid_to_idx_map_type& lid_to_elgpgp);
protected:

  // Grid name
  const std::string         m_grid_name;

  // Counters
  int                       m_num_local_elem;
  int                       m_num_gp;
  int                       m_num_vl;
  int                       m_num_local_dofs;

  // Dofs info variables
  dofs_list_type            m_dofs_gids;
  lid_to_idx_map_type       m_lid_to_idx;
};

inline FieldLayout
SEGrid::get_native_dof_layout () const
{
  using namespace ShortFieldTagsNames;

  return FieldLayout({EL,GP,GP},{m_num_local_elem,m_num_gp,m_num_gp});
}


} // namespace scream

#endif // SCREAM_SE_GRID_HPP
