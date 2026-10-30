#ifndef EAMXX_UNIQUE_GRID_REMAPPER_HPP
#define EAMXX_UNIQUE_GRID_REMAPPER_HPP

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
 * The forward direction (grid->unique grid) is a local gather: each column of the unique
 * grid must have a replica among the dofs of the source grid on the same rank.
 *
 * The backward direction (unique grid->grid) sets each dof of the source grid to the value
 * of the unique grid column with the same gid. The column may be owned by another rank, so
 * this direction involves p2p communication. The comm schedule is created at the first
 * call of remap_bwd (it's a collective operation).
 * Every dof of the source grid must correspond to a column of the unique grid (on some rank).
 *
 * Fields on the source grid can be subfields of other fields (but not subfields of
 * subfields), possibly with a dynamic subview info.
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
  void remap_bwd_impl () override;

  void initialize_device_variables ();

  // unique->grid requires communication. These methods set up/perform it.
  void setup_bwd ();
  void pack_and_send ();
  void recv_and_unpack ();

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

  // ----- unique->grid data structures ----- //

  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;

  template<typename T>
  using view_1d = typename KT::template view_ND<T,1>;
  template<typename T>
  using hview_1d = typename view_1d<T>::host_mirror_type;

  // If MpiOnDev=true, we can pass device pointers to MPI. Otherwise, we need host mirrors.
  template<typename T>
  using mpi_view_1d = std::conditional_t<MpiOnDev,view_1d<T>,hview_1d<T>>;

  bool m_bwd_setup_done = false;

  // For each unique local gid (the "overlap" grid), the list of local source dofs with such gid.
  // The dofs for the u-th unique gid are m_u2src_dofs[m_u2src_offsets[u]...m_u2src_offsets[u+1]-1]
  view_1d<int>  m_u2src_offsets;
  view_1d<int>  m_u2src_dofs;

  // The imports/exports tell us which pid we get/send unique grid columns from/to.
  // Imports are indexed by the unique local gid. Exports by the local unique grid column.
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

#endif // EAMXX_UNIQUE_GRID_REMAPPER_HPP
