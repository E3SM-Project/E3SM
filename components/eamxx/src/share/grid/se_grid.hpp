#ifndef SCREAM_SE_GRID_HPP
#define SCREAM_SE_GRID_HPP

#include "share/grid/abstract_grid.hpp"

namespace scream
{

class SEGrid : public AbstractGrid
{
public:

  // Constructor
  SEGrid (const std::string& grid_name,
          const int num_my_elements,
          const int num_gauss_pts,
          const int num_vertical_levels,
          const ekat::Comm& comm);

  virtual ~SEGrid () = default;

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  FieldLayout get_2d_scalar_layout () const override;
  FieldLayout get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const override;
  FieldLayout get_3d_scalar_layout (const bool midpoints) const override;
  FieldLayout get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const override;

  FieldTag get_partitioned_dim_tag () const override {
    return FieldTag::Element;
  }
  int get_partitioned_dim_local_size  () const override {
    return m_num_local_elem;
  }
  int get_partitioned_dim_global_size () const override {
    return m_num_global_elem;
  }

  // Retrieve list of the CG grid dofs. Const version returns a read-only field
  Field get_cg_dofs_gids ();
  Field get_cg_dofs_gids () const;

  std::shared_ptr<AbstractGrid> clone (const std::string& clone_name,
                                       const bool shallow) const override;

  bool check_valid_dofs()        const override;
  bool check_valid_lid_to_idx () const override;

protected:

  // SE dims
  int       m_num_local_elem;
  int       m_num_global_elem;
  int       m_num_gp;

  // The dofs gids for a CG version of this grid
  Field m_cg_dofs_gids;
};

} // namespace scream

#endif // SCREAM_SE_GRID_HPP
