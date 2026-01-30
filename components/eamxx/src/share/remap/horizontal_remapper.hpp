#ifndef SCREAM_HORIZONTAL_REMAPPER_HPP
#define SCREAM_HORIZONTAL_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/remap/horiz_interp_remapper_data.hpp"

#include <mpi.h>

namespace scream
{

class GridImportExport;

/*
 * A remapper to interpolate fields in the horizontal direction
 *
 * This remapper loads an interpolation sparse matrix from a map file,
 * and performs an interpolation form a src to tgt grid by means
 * of a mat-vec product. The sparse matrix encodes the interpolation
 * weights. So far, the map file is *assumed* to store the matrix in
 * triplet format, with row/col indices starting from 1.
 *
 * The remapper takes a grid and the name of the map file. The grid could
 * be either the source or the target in the remap operation, and can be
 * either the coarser or finer one.
 *
 * The mat-vec is performed in two stages
 *   1. a local mat-vec multiplication (on device)
 *   2. a pack-send-recv-unpack sequence to share data across ranks
 *
 * The order in which the two stages are performed depends on which one
 * of the two grids is coarser. The goal is to minimize the amount of data
 * that is communicated, and therefore we do MPI on the coarse grid side.
 * Therefore,
 *  - for refining remappers, we first export entries of the src field(s)
 *    that are needed on other ranks in order to perform mat-vec, and
 *    then perform the mat-vec y=W*x_overlap
 *  - for coarsening remappers, we first perform the mat-vec y_overlap=W*x
 *    and then have each rank gathering all local contributions for the
 *    entries of the tgt field it owns
 *
 * The class has to create temporaries for the intermediate fields.
 * An obvious future development would be to use some scratch memory
 * for these fields, so to not increase memory pressure.
 */

class HorizontalRemapper : public AbstractRemapper
{
public:

  HorizontalRemapper (const grid_ptr_type& src_grid,
                      const std::string& map_file,
                      const bool track_mask = false);

  ~HorizontalRemapper ();

protected:

  void registration_ends_impl () override;

  void remap_fwd_impl () override;

  // This class uses itself to remap src grid geo data to the tgt grid. But in order
  // to not pollute the remapper for later use, we must be able to clean it up after
  // remapping all the geo data.
  void clean_up ();

protected:

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void local_mat_vec (const Field& f_src, const Field& f_tgt) const;
  template<int N>
  void local_mat_vec (const Field& f_src, const Field& f_tgt, const Field& mask) const;
  template<int N>
  void rescale_masked_fields (const Field& f_tgt, const Field& f_mask) const;
  void pack_and_send ();
  void recv_and_unpack ();

protected:

  void create_ov_fields ();
  void setup_mpi_data_structures ();

  // Whether we are tracking mask fields
  bool                m_track_mask;

  // Indermediate version of the fields, on the overlap grid
  std::vector<Field>  m_ov_fields;

  // Whether each field needs to be remapped (i.e., has COL tag)
  std::vector<int>    m_needs_remap;

  // We need to keep this (and not just its content) so that the weak_ptr in HorizRemapperDataRepo
  // does not expire. This allows other remappers that need the same data to reuse it rather than
  // have the repo re-create it anew.
  std::shared_ptr<const HorizRemapperData> m_remap_data;

  // ------- MPI-related data structures -------- //

  // Offset of each field when we splice together one col of each.
  std::vector<int> m_field_offset;

  // ImportData/export info
  std::shared_ptr<GridImportExport>  m_imp_exp;

  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;

  // If MpiOnDev=true, we can pass device pointers to MPI. Otherwise, we need host mirrors.
  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;
  template<typename T>
  using hview_1d = typename view_1d<T>::HostMirror;

  template<typename T>
  using view_2d = typename KT::template view_2d<T>;

  template<typename T>
  using mpi_view_1d = std::conditional_t<MpiOnDev,view_1d<T>,hview_1d<T>>;

  // The send/recv buffers for pack/unpack operations
  view_1d<Real> m_send_buffer;
  view_1d<Real> m_recv_buffer;

  // The send/recv buf to feed to MPI (alis the above two if MpiOnDev=true)
  mpi_view_1d<Real> m_mpi_send_buffer;
  mpi_view_1d<Real> m_mpi_recv_buffer;

  // Offset of each pid's data in send/recv buffers
  view_1d<int> m_pids_send_offsets;
  view_1d<int> m_pids_recv_offsets;

  // Send/recv persistent requests
  std::vector<MPI_Request>  m_send_req;
  std::vector<MPI_Request>  m_recv_req;
};

} // namespace scream

#endif // SCREAM_HORIZONTAL_REMAPPER_HPP
