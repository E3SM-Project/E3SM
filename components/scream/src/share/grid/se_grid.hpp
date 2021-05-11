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
          const int num_global_elements,
          const int num_my_elements,
          const int num_gauss_pts,
          const int num_vertical_levels);

  virtual ~SEGrid () = default;

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  FieldLayout get_2d_scalar_layout () const override;
  FieldLayout get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const override;
  FieldLayout get_3d_scalar_layout (const bool midpoints) const override;
  FieldLayout get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const override;

  // Set the dofs gids as well as their coordinates (ie,igp,jgp) in the grid
  void set_dofs (const dofs_list_type&      dofs,
                 const lid_to_idx_map_type& lid_to_elgpgp);

  void set_geometry_data (const std::string& name, const geo_view_type& data) override;

protected:

  // SE dims
  int                       m_num_local_elem;
  int                       m_num_gp;
};

} // namespace scream

#endif // SCREAM_SE_GRID_HPP
