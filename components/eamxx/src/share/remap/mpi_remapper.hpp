#ifndef EAMXX_MPI_REMAPPER_HPP
#define EAMXX_MPI_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"
#include "share/core/eamxx_config.hpp"
#include "share/grid/grid_import_export.hpp"

#include <ekat_comm.hpp>

#include <mpi.h>
#include <memory>
#include <vector>

namespace scream
{

/*
 * Base class for remappers that exchange columns of data across ranks with p2p MPI
 *
 * It owns the send/recv buffers (and their host mirrors, if MPI cannot use device pointers)
 * and the persistent MPI requests. Derived classes decide how data is packed into/unpacked
 * from the buffers. The buffers are split in one segment per pid (sorted by pid), and each
 * segment holds a whole number of columns. A column is the splice of the columns of all the
 * registered fields (see field_col_offset). Which columns go to/come from which pid is given
 * by a GridImportExport (between a unique grid and an overlap grid).
 *
 * A typical remap does
 *   start_recvs (); pack ...; copy_send_to_mpi (); start_sends ();
 *   wait_recvs (); copy_recv_from_mpi (); unpack ...; wait_sends ();
 */
class MpiRemapper : public AbstractRemapper
{
public:
  MpiRemapper () = default;
  MpiRemapper (const grid_ptr_type& src_grid,
               const grid_ptr_type& tgt_grid)
   : AbstractRemapper(src_grid,tgt_grid) {}

protected:

  static constexpr bool MpiOnDev = SCREAM_MPI_ON_DEVICE;

  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_ND<T,1>;
  template<typename T>
  using hview_1d = typename view_1d<T>::host_mirror_type;

  // If MpiOnDev=true, we can pass device pointers to MPI. Otherwise, we need host mirrors.
  template<typename T>
  using mpi_view_1d = std::conditional_t<MpiOnDev,view_1d<T>,hview_1d<T>>;

  // Which way the data flows: from the columns of the unique grid to those of the overlap grid (we send the
  // exports of the GridImportExport, and receive the imports), or the other way around.
  enum class CommDirection { UniqueToOverlap, OverlapToUnique };

  // Allocate the buffers, and create the persistent requests.
  // If skip_self=true, no requests are created for this rank (the derived class must copy the data).
  // NOTE: this is a collective operation only in the sense that pairs of ranks must agree on the sizes.
  // This must be called after the fields are registered.
  void setup_mpi (const std::shared_ptr<const GridImportExport>& imp_exp,
                  const CommDirection dir, const bool skip_self = false);

  // The size of a column of each field (0 if the field has no COL dim), and of all fields spliced together.
  // The columns of the fields are one after the other, so the offset of a field is the sum of the sizes
  // of the previous ones.
  int field_col_size (const int i) const { return m_field_col_offset[i+1]-m_field_col_offset[i]; }
  int field_col_offset (const int i) const { return m_field_col_offset[i]; }
  int total_col_size () const { return m_field_col_offset.empty() ? 0 : m_field_col_offset.back(); }

  // The pid and local id of each column we send/receive (sorted by pid), and where each pid starts
  // (size nranks+1). These are the exports/imports of the GridImportExport, or the opposite (see CommDirection).
  using int_view_d = typename GridImportExport::view_d;
  using int_view_h = typename GridImportExport::view_h;
  int_view_d send_pids () const { return m_send_is_export ? m_imp_exp->export_pids() : m_imp_exp->import_pids(); }
  int_view_d send_lids () const { return m_send_is_export ? m_imp_exp->export_lids() : m_imp_exp->import_lids(); }
  int_view_d recv_pids () const { return m_send_is_export ? m_imp_exp->import_pids() : m_imp_exp->export_pids(); }
  int_view_d recv_lids () const { return m_send_is_export ? m_imp_exp->import_lids() : m_imp_exp->export_lids(); }
  int_view_d send_offsets () const { return m_send_is_export ? m_imp_exp->export_pid_offsets() : m_imp_exp->import_pid_offsets(); }
  int_view_d recv_offsets () const { return m_send_is_export ? m_imp_exp->import_pid_offsets() : m_imp_exp->export_pid_offsets(); }
  int_view_h send_offsets_h () const { return m_send_is_export ? m_imp_exp->export_pid_offsets_h() : m_imp_exp->import_pid_offsets_h(); }
  int_view_h recv_offsets_h () const { return m_send_is_export ? m_imp_exp->import_pid_offsets_h() : m_imp_exp->export_pid_offsets_h(); }
  int num_send_cols () const { return m_send_is_export ? m_imp_exp->num_exports() : m_imp_exp->num_imports(); }
  int num_recv_cols () const { return m_send_is_export ? m_imp_exp->num_imports() : m_imp_exp->num_exports(); }

  // Free the requests and the buffers
  void clean_mpi ();

  void start_sends ();
  void start_recvs ();
  void wait_sends ();
  void wait_recvs ();

  // No-ops if MpiOnDev=true
  void copy_send_to_mpi ();
  void copy_recv_from_mpi ();
  // Only the range [beg,end) of the recv buffer
  void copy_recv_from_mpi (const size_t beg, const size_t end);

  // The send/recv buffers for pack/unpack operations
  view_1d<Real> m_send_buffer;
  view_1d<Real> m_recv_buffer;

  // The send/recv buf to feed to MPI (alias the above two if MpiOnDev=true)
  mpi_view_1d<Real> m_mpi_send_buffer;
  mpi_view_1d<Real> m_mpi_recv_buffer;

private:

  std::shared_ptr<const GridImportExport> m_imp_exp;
  std::vector<int> m_field_col_offset;
  bool m_send_is_export = true;

  // Send/recv persistent requests
  // NOTE: remappers are used as functors in Kokkos kernels, so they get copied (and destroyed) around.
  //       Hence, the requests cannot be freed in the destructor of this class. Instead, we store them in
  //       a struct that is shared among copies, and that frees the requests when the last copy goes away.
  struct MpiRequests {
    ~MpiRequests ();
    std::vector<MPI_Request>  send;
    std::vector<MPI_Request>  recv;
  };
  std::shared_ptr<MpiRequests>  m_requests;

  void start (std::vector<MPI_Request>& reqs, const std::string& what) const;
  void wait  (std::vector<MPI_Request>& reqs, const std::string& what) const;
};

} // namespace scream

#endif // EAMXX_MPI_REMAPPER_HPP
