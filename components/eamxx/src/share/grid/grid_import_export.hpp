#ifndef EAMXX_GRID_IMPORT_EXPORT_HPP
#define EAMXX_GRID_IMPORT_EXPORT_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"       // For KokkosTypes
#include "share/util/scream_utils.hpp"  // For check_mpi_call

#include <ekat/mpi/ekat_comm.hpp>
#include <mpi.h> // We do some direct MPI calls
#include <memory>
#include <map>
#include <vector>

namespace scream
{

/*
 * Import/Export data is used to figure out where data
 * can be retrieved or sent (both in terms of remote
 * rank and remote local id) when transferring between
 * two grids. The terms import/export do not necessarily
 * imply that we are getting/sending data. Rather, they
 * have to do with which grid we're initiating the
 * transfer from. Namely,
 *   - import: the grid we have data on is potentially
 *             non unique, and the target grid is unique
 *   - export: the grid we have data on is unique, and
 *             the target grid is potentially non unique
 * Import/Export data plans can be used both for scattering
 * and gathering data. The user can use the gather/scatter
 * methods for this, but pay attention to their limitations:
 *   - they create send/recv requests at every call (no persistent requests)
 *   - they assume same data type on origin/target ranks
 *   - they operate on a particular input/output data ortanization,
 *     namely, data is organized as a map lid->vector<T>
 *   - the above point implies data must be on Host
 * These limitations imply that gather/scatter methods are only
 * for ease of use in non-performance critical code.
 * On the other hand, the import/export data (pids/lids) can
 * be used both on host and device, for more efficient pack/unpack methods.
 */

class GridImportExport {
public:
  using KT = KokkosTypes<DefaultDevice>;
  template<typename T>
  using view_1d = typename KT::view_1d<int>;

  GridImportExport (const std::shared_ptr<const AbstractGrid>& unique,
                    const std::shared_ptr<const AbstractGrid>& overlapped);
  ~GridImportExport () = default;

  template<typename T>
  void scatter (const MPI_Datatype mpi_data_t,
                const std::map<int,std::vector<T>>& src,
                      std::map<int,std::vector<T>>& dst) const;

  template<typename T>
  void gather (const MPI_Datatype mpi_data_t,
               const std::map<int,std::vector<T>>& src,
                     std::map<int,std::vector<T>>& dst) const;

  view_1d<int> num_exports_per_pid () const { return m_num_exports_per_pid; }
  view_1d<int> num_imports_per_pid () const { return m_num_imports_per_pid; }

  view_1d<int>::HostMirror num_exports_per_pid_h () const { return m_num_exports_per_pid_h; }
  view_1d<int>::HostMirror num_imports_per_pid_h () const { return m_num_imports_per_pid_h; }

  view_1d<int> import_pids () const { return m_import_pids; }
  view_1d<int> import_lids () const { return m_import_lids; }
  view_1d<int> export_pids () const { return m_export_pids; }
  view_1d<int> export_lids () const { return m_export_lids; }

  view_1d<int>::HostMirror import_pids_h () const { return m_import_pids_h; }
  view_1d<int>::HostMirror import_lids_h () const { return m_import_lids_h; }
  view_1d<int>::HostMirror export_pids_h () const { return m_export_pids_h; }
  view_1d<int>::HostMirror export_lids_h () const { return m_export_lids_h; }

protected:

  std::shared_ptr<const AbstractGrid>   m_unique;
  std::shared_ptr<const AbstractGrid>   m_overlapped;

  // All these arrays are sorted by pid. That is, all imports
  // for pid 1 come before imports for pid 2, and same for exports.
  view_1d<int>  m_import_pids;
  view_1d<int>  m_import_lids;
  view_1d<int>  m_export_pids;
  view_1d<int>  m_export_lids;

  view_1d<int>::HostMirror  m_import_pids_h;
  view_1d<int>::HostMirror  m_import_lids_h;
  view_1d<int>::HostMirror  m_export_pids_h;
  view_1d<int>::HostMirror  m_export_lids_h;

  view_1d<int>  m_num_imports_per_pid;
  view_1d<int>  m_num_exports_per_pid;

  view_1d<int>::HostMirror  m_num_imports_per_pid_h;
  view_1d<int>::HostMirror  m_num_exports_per_pid_h;

  ekat::Comm    m_comm;
};

// --------------------- IMPLEMENTATION ------------------------ //

template<typename T>
void GridImportExport::
scatter (const MPI_Datatype mpi_data_t,
         const std::map<int,std::vector<T>>& src,
               std::map<int,std::vector<T>>& dst) const
{
  using gid_type = AbstractGrid::gid_type;

  std::vector<MPI_Request> send_req, recv_req;
  
  const int nexp = m_export_lids.size();
  const int nimp = m_import_lids.size();
  auto mpi_comm = m_comm.mpi_comm();

  auto unique_gids_h = m_unique->get_dofs_gids().get_view<const gid_type*,Host>();
  auto overlap_gids_h = m_overlapped->get_dofs_gids().get_view<const gid_type*,Host>();

  // 1. Communicate to the recv pids how many items per lid
  // we need to send
  std::vector<int> send_count(m_export_lids_h.size(),0);
  for (int i=0; i<nexp; ++i) {
    auto lid = m_export_lids_h[i];
    auto gid = unique_gids_h[lid];
    send_count[i] = src.at(lid).size();
    auto& req = send_req.emplace_back();

    check_mpi_call(MPI_Isend (&send_count[i],1,MPI_INT,
                              m_export_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::scatter, creating send request (step 1)");
  }
  std::vector<int> recv_count(m_import_lids_h.size(),0);
  for (int i=0; i<nimp; ++i) {
    auto lid = m_import_lids_h[i];
    auto gid = overlap_gids_h[lid];
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv (&recv_count[i],1,MPI_INT,
                              m_import_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::scatter, creating recv request (step 1)");
  }

  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on send requests (step 1)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on recv requests (step 1)");
  send_req.clear();
  recv_req.clear();

  // 2. Size dst and create recv/send requests
  for (int i=0; i<nimp; ++i) {
    auto lid = m_import_lids_h[i];
    auto gid = overlap_gids_h[lid];
    dst[lid].resize(recv_count[i]);
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv(dst[lid].data(),recv_count[i],mpi_data_t,
                             m_import_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::scatter, creating recv request (step 2)");
  }

  for (int i=0; i<nexp; ++i) {
    auto lid = m_export_lids_h[i];
    auto gid = unique_gids_h[lid];
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend(src.at(lid).data(),src.at(lid).size(),mpi_data_t,
                             m_export_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::scatter, creating recv request (step 2)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on send requests (step 2)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on recv requests (step 2)");
  send_req.clear();
  recv_req.clear();
}

template<typename T>
void GridImportExport::
gather (const MPI_Datatype mpi_data_t,
        const std::map<int,std::vector<T>>& src,
              std::map<int,std::vector<T>>& dst) const
{
  using gid_type = AbstractGrid::gid_type;

  std::vector<MPI_Request> send_req, recv_req;
  
  const int nexp = m_export_lids.size();
  const int nimp = m_import_lids.size();
  auto mpi_comm = m_comm.mpi_comm();

  auto unique_gids_h = m_unique->get_dofs_gids().get_view<const gid_type*,Host>();
  auto overlap_gids_h = m_overlapped->get_dofs_gids().get_view<const gid_type*,Host>();

  const int num_ov_gids = overlap_gids_h.size();

  // 1. Communicate to the recv pids how many items per lid
  // we need to send
  std::vector<int> send_count(num_ov_gids,0);
  for (int i=0; i<nimp; ++i) {
    auto lid = m_import_lids_h[i];
    auto gid = overlap_gids_h[lid];
    send_count[lid] = src.at(lid).size();
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend (&send_count[lid],1,MPI_INT,
                              m_import_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::gather, creating send request (step 1)");
  }
  std::vector<int> recv_count(m_export_lids_h.size(),0);
  for (int i=0; i<nexp; ++i) {
    auto& req = recv_req.emplace_back();
    auto lid = m_export_lids_h[i];
    auto gid = unique_gids_h[lid];
    check_mpi_call(MPI_Irecv (&recv_count[i],1,MPI_INT,
                              m_export_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::gather, creating recv request (step 1)");
  }

  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on send requests (step 1)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on recv requests (step 1)");
  send_req.clear();
  recv_req.clear();

  // 2. Create recv buf and recv/send requests
  std::map<int,std::vector<T>> recv_buf;
  for (int i=0; i<nexp; ++i) {
    auto lid = m_export_lids_h[i];
    auto gid = unique_gids_h[lid];
    recv_buf[i].resize(recv_count[i]);
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv(recv_buf[i].data(),recv_count[i],mpi_data_t,
                             m_export_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::gather, creating recv request (step 2)");
  }

  for (int i=0; i<nimp; ++i) {
    auto lid = m_import_lids_h[i];
    auto& req = send_req.emplace_back();
    auto gid = overlap_gids_h[lid];
    check_mpi_call(MPI_Isend(src.at(lid).data(),src.at(lid).size(),mpi_data_t,
                             m_import_pids_h[i],gid,mpi_comm,&req),
                   "GridImportExport::gather, creating recv request (step 2)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on send requests (step 2)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on recv requests (step 2)");
  send_req.clear();
  recv_req.clear();

  // 3. Fill dst map with content from the recv buffer
  for (const auto& it : recv_buf) {
    auto lid = m_export_lids_h[it.first];
    auto& curr = dst[lid];
    auto& add = it.second;
    curr.reserve(curr.size()+add.size());
    std::move(add.begin(),add.end(),std::back_inserter(curr));
  }
}

} // namespace scream

#endif // EAMXX_GRID_IMPORT_EXPORT_DATA_HPP
