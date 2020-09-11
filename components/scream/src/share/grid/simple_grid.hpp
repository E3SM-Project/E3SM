#ifndef SCREAM_SIMPLE_HPP
#define SCREAM_SIMPLE_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"

#include "ekat/mpi/ekat_comm.hpp"

namespace scream
{

/*
 * A simple grid implementation
 *
 * This class only stores a bunch of gids, with no other additional
 * information (unlike SEGrid, which also stores info on the SE
 * element each dof belongs to.
 * This class is meant to be used mostly for testing single physics
 * parametrizations.
 */

class SimpleGrid : public AbstractGrid
{
public:

  SimpleGrid (const std::string& grid_name,
                const int num_global_columns,
                const int num_vertical_lev,
                const ekat::Comm& comm);
  virtual ~SimpleGrid () = default;

  // Grid description utilities
  GridType type () const override { return GridType::MeshFree; }
  const std::string& name () const override { return m_grid_name; }

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a 2d structured grid, this could be a set of 2 indices.
  FieldLayout get_native_dof_layout () const override;

  int get_num_vertical_levels () const override { return m_num_vl; }

  // Dofs gids utilities
  int get_num_local_dofs () const override { return m_num_my_cols; }
  const dofs_list_type& get_dofs_gids () const override { return m_dofs_gids; }
  lid_to_idx_map_type get_lid_to_idx_map () const override { return m_lid_to_idx; }

protected:

  const std::string     m_grid_name;

  int                   m_num_my_cols;
  int                   m_num_vl;
  dofs_list_type        m_dofs_gids;
  lid_to_idx_map_type   m_lid_to_idx;
};

} // namespace scream

#endif // SCREAM_SIMPLE_HPP
