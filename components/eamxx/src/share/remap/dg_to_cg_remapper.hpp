#ifndef EAMXX_DG_TO_CG_REMAPPER_HPP
#define EAMXX_DG_TO_CG_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/remap/field_views_repo.hpp"
#include "share/core/eamxx_config.hpp"

#include "share/grid/grid_import_export.hpp"

#include <ekat_pack.hpp>

#include <mpi.h>
#include <memory>
#include <vector>

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
 * The forward direction (DG->CG) is a local gather: each CG column must correspond
 * to a dof on the same rank.
 *
 * The backward direction (CG->DG) sets each dof of the SE grid to the value of the
 * CG column with the same gid. The CG column may be owned by another rank, so this
 * direction involves p2p communication. The comm schedule is created at the first
 * call of remap_bwd (it's a collective operation).
 * Every dof of the SE grid must correspond to a column of the CG grid (on some rank).
 *
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
  void remap_bwd_impl () override;

  void initialize_device_variables ();

  // CG->DG requires communication. These methods set up/perform it.
  void setup_bwd ();
  void pack_and_send ();
  void recv_and_unpack ();

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

  // ----- CG->DG data structures ----- //

  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;

  template<typename T>
  using view_1d = typename KT::template view_ND<T,1>;
  template<typename T>
  using hview_1d = typename view_1d<T>::host_mirror_type;

  // If MpiOnDev=true, we can pass device pointers to MPI. Otherwise, we need host mirrors.
  template<typename T>
  using mpi_view_1d = std::conditional_t<MpiOnDev,view_1d<T>,hview_1d<T>>;

  bool m_bwd_setup_done = false;

  // For each unique local CG gid (the "overlap" grid), the list of local SE dofs with such gid.
  // The dofs for the u-th CG gid are m_u2se_dofs[m_u2se_offsets[u]...m_u2se_offsets[u+1]-1]
  view_1d<int>  m_u2se_offsets;
  view_1d<int>  m_u2se_dofs;

  // The imports/exports tell us which pid we get/send CG columns from/to.
  // Imports are indexed by the unique local CG gid. Exports by the local CG column.
  // All arrays are sorted by pid, and the offsets (size nranks+1) tell us where each pid starts.
  int m_num_imports = 0;
  int m_num_exports = 0;
  view_1d<int>  m_imp_pids;
  view_1d<int>  m_imp_lids;
  view_1d<int>  m_imp_offsets;
  view_1d<int>  m_exp_pids;
  view_1d<int>  m_exp_lids;
  view_1d<int>  m_exp_offsets;

  // Size of one column of each field (and of all fields together), and offset of each field in the latter.
  int m_total_col_size = 0;
  view_1d<int>  m_f_col_size;
  view_1d<int>  m_f_col_offset;

  // The send/recv buffers for pack/unpack operations
  view_1d<Real> m_send_buffer;
  view_1d<Real> m_recv_buffer;

  // The send/recv buf to feed to MPI (alias the above two if MpiOnDev=true)
  mpi_view_1d<Real> m_mpi_send_buffer;
  mpi_view_1d<Real> m_mpi_recv_buffer;

  // Send/recv persistent requests
  // NOTE: this class is used as a functor in Kokkos kernels, so it gets copied (and destroyed) around.
  //       Hence, the requests cannot be freed in the destructor of this class. Instead, we store them in
  //       a struct that is shared among copies, and that frees the requests when the last copy goes away.
  struct MpiRequests {
    ~MpiRequests ();
    std::vector<MPI_Request>  send;
    std::vector<MPI_Request>  recv;
  };
  std::shared_ptr<MpiRequests>  m_requests;

  template <typename MT>
  KOKKOS_FUNCTION
  void local_remap_2d (const MT& team) const;

  template <typename ScalarT, typename MT>
  KOKKOS_FUNCTION
  void local_remap_3d (const MT& team) const;

public:
  // These should be morally private, but CUDA complains that they cannot be private/protected
  struct RemapTag {};
  struct PackTag {};
  struct UnpackTag {};

  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const RemapTag&, const MT &team) const;
  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const PackTag&, const MT &team) const;
  template<typename MT>
  KOKKOS_INLINE_FUNCTION
  void operator()(const UnpackTag&, const MT &team) const;
};

} // namespace scream

#endif // EAMXX_DG_TO_CG_REMAPPER_HPP
