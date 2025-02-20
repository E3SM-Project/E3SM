#ifndef SCREAM_REFINING_REMAPPER_RMA_HPP
#define SCREAM_REFINING_REMAPPER_RMA_HPP

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/grid/remap/horiz_interp_remapper_base.hpp"
#include "share/util/eamxx_utils.hpp"
#include "eamxx_config.h"

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

class RefiningRemapperRMA : public HorizInterpRemapperBase
{
public:

  RefiningRemapperRMA (const grid_ptr_type& tgt_grid,
                       const std::string& map_file);

  ~RefiningRemapperRMA ();

protected:

  void remap_fwd_impl () override;

  void setup_mpi_data_structures () override;

  // This class uses itself to remap src grid geo data to the tgt grid. But in order
  // to not pollute the remapper for later use, we must be able to clean it up after
  // remapping all the geo data.
  void clean_up ();

  // Wrap a pointer in an MPI_Win
  template<typename T>
  MPI_Win get_mpi_window (T* v, int n) const {
    MPI_Win win;
    check_mpi_call (MPI_Win_create(v,n*sizeof(T),sizeof(T),
                                   MPI_INFO_NULL,m_comm.mpi_comm(),&win),
                    "MPI_Win_create");
    return win;
  }

  MPI_Group             m_mpi_group = MPI_GROUP_NULL;

  // Unfortunately there is no GPU-aware mpi for RMA operations.
  //static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;
  static constexpr bool MpiOnDev = false;

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
