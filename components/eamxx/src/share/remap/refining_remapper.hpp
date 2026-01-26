#ifndef SCREAM_REFINING_REMAPPER_HPP
#define SCREAM_REFINING_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/remap/horiz_interp_remapper_base.hpp"

#include <ekat_pack.hpp>

#include <mpi.h>

namespace scream
{

class GridImportExport;

/*
 * A remapper to interpolate fields on a coarser grid
 *
 * This remapper loads an interpolation sparse matrix from a map file,
 * and performs an interpolation form a coarse to a fine grid by means
 * of a mat-vec product. The sparse matrix encodes the interpolation
 * weights. So far, the map file is *assumed* to store the matrix in
 * triplet format, with row/col indices starting from 1.
 *
 * The remapper takes a src grid and the name of the map file. From here,
 * it creates the tgt grid, and all the internal structures needed for
 * an efficient mat-vec product at runtime.
 *
 * The mat-vec is performed in two stages:
 *   1. Perform a pack-send-recv-unpack sequence via MPI, to ensure
 *      all ranks have access to the src field(s) entries they need
 *      for the mat-vec.
 *   2. Perform a local mat-vec multiplication (on device)
 *
 * The class has to create temporaries for the intermediate fields.
 * An obvious future development would be to use some scratch memory
 * for these fields, so to not increase memory pressure.
 *
 */

class RefiningRemapper : public HorizInterpRemapperBase
{
public:

  RefiningRemapper (const grid_ptr_type& tgt_grid,
                    const std::string& map_file);

  ~RefiningRemapper ();

protected:

  void remap_fwd_impl () override;

  void setup_mpi_data_structures () override;

  void clean_up ();

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void pack_and_send ();
  void recv_and_unpack ();

protected:

  // If MpiOnDev=true, we pass device pointers to MPI. Otherwise, we use host mirrors.
  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;
  template<typename T>
  using mpi_view_1d = typename std::conditional<
                        MpiOnDev,
                        view_1d<T>,
                        typename view_1d<T>::HostMirror
                      >::type;

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

#endif // SCREAM_REFINING_REMAPPER_HPP
