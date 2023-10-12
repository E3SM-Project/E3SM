#ifndef SCREAM_REFINING_REMAPPER_RMA_HPP
#define SCREAM_REFINING_REMAPPER_RMA_HPP

#include "share/grid/remap/abstract_remapper.hpp"
#include "scream_config.h"

#include "ekat/ekat_pack.hpp"

#include <mpi.h>

namespace scream
{

/*
 * A remapper to interpolate fields on a finer grid
 *
 * This remapper loads an interpolation sparse matrix from a map file,
 * and performs an interpolation form a coarse to a fine grid by means
 * of a mat-vec product. The sparse matrix encodes the interpolation
 * weights. So far, the map file is *assumed* to store the matrix in
 * triplet format, with row/col indices starting from 1.
 *
 * The remapper takes a tgt grid and the name of the map file. From here,
 * it creates the src grid, and all the internal structures needed for
 * an efficient mat-vec product at runtime.
 *
 * The mat-vec is performed in two stages:
 *   1. Import remote entries of the source fields into an overlapped
 *      partition, so that each rank has all the entries it needs to
 *      performe a local mat-vec product.
 *   2. Perform the local mat-vec product (on device), producing using
 *      as input fields the onese produced by step one.
 *
 * The class has to create temporaries for the intermediate fields.
 * An obvious future development would be to use some scratch memory
 * for these fields, so to not increase memory pressure.
 *
 * All the MPI operations performed by this class are implemented with
 * one-sided (or RMA) MPI routines. One-sided MPI has been in the MPI
 * standard since 2.0, but its support is still sub-optimal, due to
 * limited effort in optimizing it by the vendors. Furthermore, as of
 * Oct 2023, RMA operations are not supported by GPU-aware implementations.
 */

class RefiningRemapperRMA : public AbstractRemapper
{
public:

  RefiningRemapperRMA (const grid_ptr_type& tgt_grid,
                    const std::string& map_file);

  ~RefiningRemapperRMA ();

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
    EKAT_ERROR_MSG ("RefiningRemapperRMA only supports fwd remapping.\n");
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

protected:
  void check_mpi_call (int err, const std::string& context) const;

  struct Triplet {
    gid_type row;
    gid_type col;
    Real  w;
  };

  std::vector<Triplet>
  get_my_triplets (const std::string& map_file,
                   const grid_ptr_type& src_grid);

  // Wrap a pointer in an MPI_Win
  template<typename T>
  MPI_Win get_mpi_window (T* v, int n) const {
    MPI_Win win;
    check_mpi_call (MPI_Win_create(v,n*sizeof(T),sizeof(T),
                                   MPI_INFO_NULL,m_comm.mpi_comm(),&win),
                    "MPI_Win_create");
    return win;
  }

  ekat::Comm            m_comm;
  MPI_Group             m_mpi_group = MPI_GROUP_NULL;

  // Unfortunately there is no GPU-aware mpi for RMA operations.
  //static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;
  static constexpr bool MpiOnDev = false;

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

  // ------- MPI data structures -------- //

  // For each GID in m_ov_src_grid, store the pid it belongs
  // to in m_src_grid, and the local id on that pid.
  std::vector<int>          m_remote_pids;
  std::vector<int>          m_remote_lids;

  // Column info for each field.
  // Notes:
  //  - for subfields, col_stride!=col_size, otherwise they match
  //  - col_offset!=0 only for subfield that are not the 0-th entry along subf dim.
  //  - in general, col_data = col_stride*icol+col_offset.
  //  - strides/offsets  are *only* for m_src_fields (ov_src are contiguous, and tgt are only
  //    accessed via get_view).
  std::vector<int>          m_col_size;
  std::vector<int>          m_col_stride;
  std::vector<int>          m_col_offset;

  // One MPI window object for each field
  std::vector<MPI_Win>      m_mpi_win;
};

} // namespace scream

#endif // SCREAM_REFINING_REMAPPER_RMA_HPP
