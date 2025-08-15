#ifndef SCREAM_COARSENING_REMAPPER_HPP
#define SCREAM_COARSENING_REMAPPER_HPP

#include "share/grid/remap/horiz_interp_remapper_base.hpp"
#include "eamxx_config.h"

#include <mpi.h>

namespace scream
{

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
 * The setup as well as the runtime operations use classic send/recv
 * MPI calls, where data is packed in a buffer and sent to the recv rank,
 * where it is then unpacked and accumulated into the result.
 */

class CoarseningRemapper : public HorizInterpRemapperBase
{
public:

  CoarseningRemapper (const grid_ptr_type& src_grid,
                      const std::string& map_file,
                      const bool track_mask = false,
                      const bool populate_tgt_grid_geo_data = true);

  ~CoarseningRemapper ();

protected:

  void registration_ends_impl () override;

  void remap_fwd_impl () override;

  template<typename T>
  using view_2d = typename KT::template view_2d<T>;

  void setup_mpi_data_structures () override;

  std::vector<int> get_pids_for_recv (const std::vector<int>& send_to_pids) const;

  std::map<int,std::vector<int>>
  recv_gids_from_pids (const std::map<int,std::vector<int>>& pid2gids_send) const;

  // This class uses itself to remap src grid geo data to the tgt grid. But in order
  // to not pollute the remapper for later use, we must be able to clean it up after
  // remapping all the geo data.
  void clean_up ();

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void local_mat_vec (const Field& f_src, const Field& f_tgt, const Field& mask) const;
  template<int N>
  void rescale_masked_fields (const Field& f_tgt, const Field& f_mask) const;
  void pack_and_send ();
  void recv_and_unpack ();
  // Overload, not hide
  using HorizInterpRemapperBase::local_mat_vec;

protected:

  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;

  // If MpiOnDev=true, we can pass device pointers to MPI. Otherwise, we need host mirrors.
  template<typename T>
  using mpi_view_1d = typename std::conditional<
                        MpiOnDev,
                        view_1d<T>,
                        typename view_1d<T>::HostMirror
                      >::type;

  // Mask fields, if needed
  bool                  m_track_mask;

  // ------- MPI data structures -------- //

  // The send/recv buf for pack/unpack
  view_1d<Real>         m_send_buffer;
  view_1d<Real>         m_recv_buffer;

  // The send/recv buf to feed to MPI.
  // If MpiOnDev=true, they simply alias the ones above
  mpi_view_1d<Real>     m_mpi_send_buffer;
  mpi_view_1d<Real>     m_mpi_recv_buffer;

  // Offset of each field on each PID in send/recv buffers.
  // E.g., offset(3,2)=10 means the offset of data from field 3
  //       to be sent to PID 2 is 10.
  view_2d<int>          m_send_f_pid_offsets;
  view_2d<int>          m_recv_f_pid_offsets;

  // Reorder the lids so that all lids to send to PID n
  // come before those for PID N+1. The meaning is
  //   lids_pids(i,0) = ith lid to send to PID=lids_pids(i,1)
  // Note: send lids are the lids of gids in the ov_tgt_grid.
  //       But here, dofs are ordered differently, so that all dofs to
  //       send to the same PID are contiguous.
  view_2d<int>          m_send_lids_pids;

  // Store the start of lids to send to each PID in the view above
  view_1d<int>          m_send_pid_lids_start;

  // Unlike the packing for sends, unpacking after the recv can cause
  // race conditions. Hence, we ||ize of tgt lids, and process separate
  // contributions from separate PIDs serially. To do so, we use the
  // pidpos array that, for each contribution, it tells us the PID where it
  // comes from, and the position of that dof in the list of dofs that
  // this PID sends us. We sort this array by the lid of the dof, and
  // use the beg/end arrays to store the begin/end of contributions for
  // each tgt lid inside the pidpos view.
  // E.g., if beg(2)=10 and end(2)=13, then tgt lid 2 has three contributions,
  // stored in rows 10,11,and 12 of the pidpos view. Furthermore, if, say,
  // pidpos(11,0)=3, and pidpos(11,1)=19, then the 2nd contribution for lid=2
  // comes from pid 3, and pid 3 packed that dof as the 19th in the list of
  // dofs it sent to us.
  view_2d<int>          m_recv_lids_pidpos;
  view_1d<int>          m_recv_lids_beg;
  view_1d<int>          m_recv_lids_end;

  // Send/recv requests
  std::vector<MPI_Request>  m_recv_req;
  std::vector<MPI_Request>  m_send_req;
};

} // namespace scream

#endif // SCREAM_COARSENING_REMAPPER_HPP
