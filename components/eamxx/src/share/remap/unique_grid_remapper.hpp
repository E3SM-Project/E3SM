#ifndef EAMXX_UNIQUE_GRID_REMAPPER_HPP
#define EAMXX_UNIQUE_GRID_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/remap/field_views_repo.hpp"
#include "share/core/eamxx_config.hpp"

#include <ekat_pack.hpp>

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
 * communication across ranks. Fields on the source grid can be subfields of other fields
 * (but not subfields of subfields), possibly with a dynamic subview info.
 */

class UniqueGridRemapper : public AbstractRemapper
{
public:
  using device_type = FieldViewsRepo::device_type;
  using KT          = FieldViewsRepo::KT;
  using pack_type   = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  // The target grid is the unique grid of the source grid (see AbstractGrid::get_unique_grid)
  explicit UniqueGridRemapper (const grid_ptr_type& grid);

  ~UniqueGridRemapper () = default;

  bool compatible_layouts (const FieldLayout& src,
                           const FieldLayout& tgt) const override {
    return src.type()==tgt.type();
  }

  bool is_valid_src_layout (const FieldLayout& layout) const override;

protected:

  void registration_ends_impl () override;
  void remap_fwd_impl () override;

  void initialize_device_variables ();

  int m_num_unique_cols;

  // For each dof in the source grid, the indices that address it in a field (K of them: 1 or 3)
  typename Field::view_dev_t<const int**>  m_src_lid2idx;
  int m_src_K;

  // For each column in the unique grid, the corresponding local dof in the source grid
  KT::view_ND<int,1>  m_unique2src;

  FieldViewsRepo  m_src_repo;
  FieldViewsRepo  m_tgt_repo;

  enum AllocPropType : int {
    PackAlloc = 0,
    RealAlloc = 2
  };

  // Info on the fields (all identical for src and tgt fields)
  KT::view_ND<int,1>  m_layout;
  KT::view_ND<int,1>  m_pack_alloc_property;

  // Only meaningful for 3d fields. Set to -1 for 2d fields, in case wrongfully used
  KT::view_ND<int,1>  m_num_levels;

  // Only meaningful for vector fields. Set to -1 for scalar fields, in case wrongfully used.
  // Note: we can't use the view extents, since the last extent may be padded
  KT::view_ND<int,1>  m_vec_dim;

  template <typename MT>
  KOKKOS_FUNCTION
  void local_remap_2d (const MT& team) const;

  template <typename ScalarT, typename MT>
  KOKKOS_FUNCTION
  void local_remap_3d (const MT& team) const;

public:
  // These should be morally private, but CUDA complains that they cannot be private/protected
  struct RemapTag {};

  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const RemapTag&, const MT &team) const;
};

} // namespace scream

#endif // EAMXX_UNIQUE_GRID_REMAPPER_HPP
