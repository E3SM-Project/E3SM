#ifndef EAMXX_DG_TO_CG_REMAPPER_HPP
#define EAMXX_DG_TO_CG_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/remap/field_views_repo.hpp"
#include "share/core/eamxx_config.hpp"

#include <ekat_pack.hpp>

namespace scream
{

// For each column in the CG grid, find the local index of a dof in the SE grid
// with the same gid. Such dof is looked up among the *unique* (CG) dofs of the SE grid.
// Notice that such dof may not be unique (if the column is on an edge of a SE element),
// but we don't care: since SE fields are continuous, any match will do.
// An error is raised if some column cannot be matched.
KokkosTypes<DefaultDevice>::view_1d<int>
compute_cg_to_se_dof_map (const AbstractGrid& se_grid,
                          const AbstractGrid& cg_grid);

/*
 * A remapper from the (DG) SE grid to its CG representation
 *
 * A SE grid stores the dofs of each element separately (DG), so dofs on the
 * edge of an element are duplicated. The target grid is a point grid whose
 * columns are unique dofs of the SE grid, that is, whose gids appear among the
 * CG gids of the SE grid (e.g., the "physics_gll" grid in Homme). Each column
 * in the target grid receives the value of the corresponding dof in the source
 * grid. No interpolation is involved.
 *
 * Only the forward direction (SE->CG) is supported, since the opposite
 * direction would require some form of communication across elements.
 * Fields on the SE grid can be subfields of other fields (but not subfields
 * of subfields), possibly with a dynamic subview info.
 */

class DGToCGRemapper : public AbstractRemapper
{
public:
  using device_type = FieldViewsRepo::device_type;
  using KT          = FieldViewsRepo::KT;
  using pack_type   = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  DGToCGRemapper (const grid_ptr_type& se_grid,
                  const grid_ptr_type& cg_grid);

  ~DGToCGRemapper () = default;

  bool compatible_layouts (const FieldLayout& src,
                           const FieldLayout& tgt) const override {
    return src.type()==tgt.type();
  }

  bool is_valid_src_layout (const FieldLayout& layout) const override;

protected:

  void registration_ends_impl () override;
  void remap_fwd_impl () override;

  void initialize_device_variables ();

  int m_num_cg_cols;

  // For each dof in the SE grid, its (element,igp,jgp) indices
  typename Field::view_dev_t<const int**>  m_lid2elgp;

  // For each column in the CG grid, the corresponding local dof in the SE grid
  KT::view_ND<int,1>  m_cg2se;

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

#endif // EAMXX_DG_TO_CG_REMAPPER_HPP
