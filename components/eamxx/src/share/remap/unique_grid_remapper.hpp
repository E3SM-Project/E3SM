#ifndef EAMXX_UNIQUE_GRID_REMAPPER_HPP
#define EAMXX_UNIQUE_GRID_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/core/eamxx_config.hpp"

#include <vector>

namespace scream
{

// For each column of the unique grid, find the local dof of the (source) grid that is a replica
// of it, that is, whose entry in the unique grid gids is the gid of the column.
// If more than one dof matches (replicated dofs), any match will do, since replicas
// hold the same value. An error is raised if some column cannot be matched.
KokkosTypes<DefaultDevice>::view_1d<int>
compute_unique_grid_dof_map (const AbstractGrid& grid,
                             const AbstractGrid& unique_grid);

/*
 * A remapper from a grid to its unique grid
 *
 * Some grids have distinct dofs that represent the same point (replicated dofs). E.g.,
 * a SE grid stores the dofs of each element separately, so dofs on the edge of an element
 * are repeated. The unique grid collapses replicated dofs (see AbstractGrid::get_unique_grid).
 * Each column in the target (unique) grid receives the value of one of the replicas
 * in the source grid that map to it. No interpolation is involved.
 *
 * The source grid is addressed through its lid->idx map, which can have 1 index per dof
 * (e.g., (col) layouts) or 3 (e.g., (elem,gp,gp) layouts). The target grid must be
 * addressed with 1 index per dof. Every column of the target grid must have a
 * replica among the source grid dofs on the same rank.
 *
 * Only the forward direction is supported, since the opposite direction would require
 * communication across ranks. Fields can be subfields of other fields, possibly with a
 * dynamic subview info (the views are extracted at every remap).
 */

class UniqueGridRemapper : public AbstractRemapper
{
public:
  using device_type = Field::device_t;
  using KT          = KokkosTypes<device_type>;

  // The target grid is the unique grid of the source grid (see AbstractGrid::get_unique_grid)
  explicit UniqueGridRemapper (const grid_ptr_type& grid);

  ~UniqueGridRemapper () = default;

  bool compatible_layouts (const FieldLayout& src,
                           const FieldLayout& tgt) const override {
    // The layouts must be congruent, except for the horizontal dimension(s), which
    // are different by design (the kernels use the source dims to loop over the data)
    if (src.type()!=tgt.type() or src.rank()==0)
      return false;
    if (src.is_vector_layout() and src.get_vector_dim()!=tgt.get_vector_dim())
      return false;
    const bool is_3d = src.type()==LayoutType::Scalar3D or src.type()==LayoutType::Vector3D;
    if (is_3d and (src.tags().back()!=tgt.tags().back() or src.dims().back()!=tgt.dims().back()))
      return false;
    return true;
  }

  bool is_valid_src_layout (const FieldLayout& layout) const override;

protected:

  void registration_ends_impl () override;
  void remap_fwd_impl () override;

  // Re-target the component subfields of the i-th field (if it is a vector field) to the c-th component
  void set_component (const int i, const int c);

  int m_num_unique_cols;

  // For each dof in the source grid, the indices that address it in a field (K of them: 1 or 3)
  typename Field::view_dev_t<const int**>  m_src_lid2idx;
  int m_src_K;

  // For each column in the unique grid, the corresponding local dof in the source grid
  KT::view_ND<int,1>  m_unique2src;

  // Vector fields are remapped one component at a time. For each registered field, we store
  // whether it is a vector field, its number of components (1 if not), and the scalar fields
  // that the kernels work with: the (dynamic) component subfields for vector fields, and the
  // fields themselves otherwise.
  std::vector<bool>   m_is_vec;
  std::vector<int>    m_ncomp;
  std::vector<Field>  m_src_scalar;
  std::vector<Field>  m_tgt_scalar;
};

} // namespace scream

#endif // EAMXX_UNIQUE_GRID_REMAPPER_HPP
