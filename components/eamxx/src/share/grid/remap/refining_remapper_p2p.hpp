#ifndef SCREAM_REFINING_REMAPPER_P2P_HPP
#define SCREAM_REFINING_REMAPPER_P2P_HPP

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/remap/horiz_interp_remapper_base.hpp"
#include "scream_config.h"

#include "ekat/ekat_pack.hpp"

#include <mpi.h>

namespace scream
{

class GridImportExport;

/*
 * A remapper to interpolate fields on a coarser grid
 *
 * This remapper loads an interpolation sparse matrix from a map file,
 * and performs an interpolation form a fine to a coarse grid by means
 * of a mat-vec product. The sparse matrix encodes the interpolation
 * weights. So far, the map file is *assumed* to store the matrix in
 * triplet format, with row/col indices starting from 1.
 *
 * The remapper takes a src grid and the name of the map file. From here,
 * it creates the tgt grid, and all the internal structures needed for
 * an efficient mat-vec product at runtime.
 *
 * The mat-vec is performed in two stages:
 *   1. Perform a local mat-vec multiplication (on device), producing intermediate
 *      output fields that have "duplicated" entries (that is, 2+ MPI
 *      ranks could all own a piece of the result for the same dof).
 *   2. Perform a pack-send-recv-unpack sequence via MPI, to accumulate
 *      partial results on the rank that owns the dof in the tgt grid.
 *
 * The class has to create temporaries for the intermediate fields.
 * An obvious future development would be to use some scratch memory
 * for these fields, so to not increase memory pressure.
 *
 * The setup of the class uses a bunch of RMA mpi operations, since they
 * are more convenient when ranks don't know where data is coming from
 * or how much data is coming from each rank. The runtime operations,
 * however, use the classic send/recv paradigm, where data is packed in
 * a buffer, sent to the recv rank, and then unpacked and accumulated
 * into the result.
 */

class RefiningRemapperP2P : public AbstractRemapper,
                            public HorizInterpRemapperBase
{
public:

  RefiningRemapperP2P (const grid_ptr_type& tgt_grid,
                    const std::string& map_file);

  ~RefiningRemapperP2P ();

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override;
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override;

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const override {
    // Same type of layout, and same sizes except for possibly the first one
    constexpr auto COL = ShortFieldTagsNames::COL;
    return get_layout_type(src.tags())==get_layout_type(tgt.tags()) &&
           src.strip_dim(COL)==tgt.strip_dim(COL);
  }

protected:

  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_src_fields[ifield].get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_tgt_fields[ifield].get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_src_fields[ifield];
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_tgt_fields[ifield];
  }

  void do_registration_begins () override { /* Nothing to do here */ }

  void do_register_field (const identifier_type& src, const identifier_type& tgt) override;

  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) override;

  void do_registration_ends () override;

  void do_remap_fwd () override;

  void do_remap_bwd () override {
    EKAT_ERROR_MSG ("RefiningRemapperP2P only supports fwd remapping.\n");
  }

protected:

  using KT = KokkosTypes<DefaultDevice>;
  using gid_type = AbstractGrid::gid_type;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;

  void create_ov_src_fields ();
  void setup_mpi_data_structures ();

  // This class uses itself to remap src grid geo data to the tgt grid. But in order
  // to not pollute the remapper for later use, we must be able to clean it up after
  // remapping all the geo data.
  void clean_up ();

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void local_mat_vec (const Field& f_src, const Field& f_tgt) const;
  void pack_and_send ();
  void recv_and_unpack ();

protected:
  // void check_mpi_call (int err, const std::string& context) const;

  ekat::Comm            m_comm;

  // If MpiOnDev=true, we pass device pointers to MPI. Otherwise, we use host mirrors.
  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;
  template<typename T>
  using mpi_view_1d = typename std::conditional<
                        MpiOnDev,
                        view_1d<T>,
                        typename view_1d<T>::HostMirror
                      >::type;

  // An "overlapped" src grid, that is a version of the src grid where
  // ranks own all cols that are affecting local dofs in their tgt grid
  grid_ptr_type         m_ov_src_grid;

  // Source, target, and overlapped-target fields
  std::vector<Field>    m_src_fields;
  std::vector<Field>    m_ov_src_fields;
  std::vector<Field>    m_tgt_fields;

  // ----- Sparse matrix CRS representation ---- //
  view_1d<int>    m_row_offsets;
  view_1d<int>    m_col_lids;
  view_1d<Real>   m_weights;

  // ----- Data structures for pack/unpack and MPI ----- //

  // Exclusive scan sum of the col size of each field
  std::vector<int> m_fields_col_sizes_scan_sum;

  // ImportData/export info
  std::shared_ptr<GridImportExport>  m_imp_exp;

  // The send/recv buffers for pack/unpack operations
  view_1d<Real> m_send_buffer;
  view_1d<Real> m_recv_buffer;

  // The send/recv buf to feed to MPI.
  // If MpiOnDev=true, they simply alias the ones above
  mpi_view_1d<Real>     m_mpi_send_buffer;
  mpi_view_1d<Real>     m_mpi_recv_buffer;

  // Offset of each pid in send/recv buffers
  view_1d<int>  m_pids_send_offsets;
  view_1d<int>  m_pids_recv_offsets;

  // For each col, its position within the set of cols
  // sent/recv to/from the corresponding remote
  view_1d<int>  m_send_col_pos;
  view_1d<int>  m_recv_col_pos;

  // Send/recv persistent requests
  std::vector<MPI_Request>  m_send_req;
  std::vector<MPI_Request>  m_recv_req;
};

} // namespace scream

#endif // SCREAM_REFINING_REMAPPER_P2P_HPP
