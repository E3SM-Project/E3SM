#ifndef EAMXX_GRID_IMPORT_EXPORT_HPP
#define EAMXX_GRID_IMPORT_EXPORT_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/eamxx_types.hpp"       // For KokkosTypes
#include "share/util/eamxx_utils.hpp"  // For check_mpi_call

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

  using gid_type = AbstractGrid::gid_type;

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

// Both gather and scatter proceed in 4 stages. If pid1 needs to send data to pid2, then
//  1. pid1 sends pid2 the list of GIDs [gid1,...,gidN] it will send. Although pid2
//     already knows the GIDs, this step specifies the ORDER, since GIDs
//     may be ordered differently on different pids
//  2. pid1 sends pid2 the number of T's it will send for each of those gids
//  3. pid1 sends pid3 all T's (first T's for gid1, then T's for gid2, etc)
//  4. Finally, pid2 parses the data received, stuffing it into dst, according
//     to the lid on pid2 (needs converting gids to lids)
// The two only differ in which members to use. They are mostly specular: where one
// uses m_unique grid, the other uses m_overlapped grid, and where one use m_export_lids/pids,
// the other uses m_import_lids/pids. That's b/c they do the same operation, just
// in opposite directions.

template<typename T>
void GridImportExport::
scatter (const MPI_Datatype mpi_data_t,
         const std::map<int,std::vector<T>>& src,
               std::map<int,std::vector<T>>& dst) const
{
  const auto tag = 0;

  std::vector<MPI_Request> send_req, recv_req;
  
  const int nexp = m_export_lids.size();
  const int nimp = m_import_lids.size();
  auto mpi_comm = m_comm.mpi_comm();
  auto mpi_gid_t = ekat::get_mpi_type<gid_type>();

  auto gids_h = m_unique->get_dofs_gids().get_view<const gid_type*,Host>();

  // 1. Communicate GIDs lists to recv pids
  std::map<int,std::vector<gid_type>> send_pid2gids;
  for (int i=0; i<nexp; ++i) {
    auto pid = m_export_pids_h[i];
    auto lid = m_export_lids_h[i];
    send_pid2gids[pid].push_back(gids_h[lid]);
  }
  for (auto& [pid,gids] : send_pid2gids) {
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend (gids.data(),gids.size(),mpi_gid_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::scatter, creating send request (step 1)");
  }
  std::map<int,std::vector<gid_type>> recv_pid2gids;
  for (int i=0; i<nimp; ++i) {
    auto pid = m_import_pids_h[i];
    recv_pid2gids[pid].push_back(-1);
  }
  for (auto& [pid,gids] : recv_pid2gids) {
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv (gids.data(),gids.size(),mpi_gid_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::scatter, creating recv request (step 1)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on send requests (step 1)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on recv requests (step 1)");
  send_req.clear();
  recv_req.clear();
  
  // 2. Communicate T's count for each GID to recv pids
  std::map<int,std::vector<int>> send_pid2count;
  for (int i=0; i<nexp; ++i) {
    auto pid = m_export_pids_h[i];
    auto lid = m_export_lids_h[i];
    send_pid2count[pid].push_back(src.at(lid).size());
  }
  for (auto& [pid,count] : send_pid2count) {
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend (count.data(),count.size(),MPI_INT,pid,tag,mpi_comm,&req),
                   "GridImportExport::scatter, creating send request (step 2)");
  }

  std::map<int,std::vector<int>>      recv_pid2count;
  for (auto& [pid,gids] : recv_pid2gids) {
    recv_pid2count[pid].resize(gids.size());
  }
  for (auto& [pid,count] : recv_pid2count) {
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv (count.data(),count.size(),MPI_INT,pid,tag,mpi_comm,&req),
                   "GridImportExport::scatter, creating recv request (step 2)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on send requests (step 2)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on recv requests (step 2)");
  send_req.clear();
  recv_req.clear();

  // 3. Pack and send the data.
  std::map<int,std::vector<T>> send_pid2data;
  const auto& gid2lid = m_unique->get_gid2lid_map();
  for (const auto& [pid,gids] : send_pid2gids) {
    auto& data = send_pid2data[pid];
    for (auto g : gids) {
      auto lid = gid2lid.at(g);
      data.insert(data.end(),src.at(lid).begin(),src.at(lid).end());
    }
  }
  for (auto& [pid,data] : send_pid2data) {
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend (data.data(),data.size(),mpi_data_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::scatter, creating send request (step 3)");
  }

  std::map<int,std::vector<T>> recv_pid2data;
  for (const auto& [pid,count] : recv_pid2count) {
    int n = 0;
    for (auto c : count) n += c;

    auto& data = recv_pid2data[pid];
    data.resize(n);
  }
  for (auto& [pid,data] : recv_pid2data) {
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv (data.data(),data.size(),mpi_data_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::scatter, creating recv request (step 3)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on send requests (step 3)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::scatter, waiting on recv requests (step 3)");
  send_req.clear();
  recv_req.clear();

  // 4. Unpack received data in dst map
  const auto& recv_gid2lid = m_overlapped->get_gid2lid_map();
  for (const auto& [pid,data] : recv_pid2data) {
    const auto& count = recv_pid2count[pid];
    const auto& gids  = recv_pid2gids[pid];
    const int num_gids = count.size();
    for (int i=0, pos=0; i<num_gids; ++i) {
      auto lid = recv_gid2lid.at(gids[i]);
      auto curr_sz = dst[lid].size();
      dst[lid].resize(curr_sz+count[i]);
      for (int k=0; k<count[i]; ++k, ++pos) {
        dst[lid][curr_sz+k] = data[pos];
      }
    }
  }
}

template<typename T>
void GridImportExport::
gather (const MPI_Datatype mpi_data_t,
        const std::map<int,std::vector<T>>& src,
              std::map<int,std::vector<T>>& dst) const
{
  const auto tag = 0;

  std::vector<MPI_Request> send_req, recv_req;
  
  const int nexp = m_export_lids.size();
  const int nimp = m_import_lids.size();
  auto mpi_comm = m_comm.mpi_comm();
  auto mpi_gid_t = ekat::get_mpi_type<gid_type>();

  auto ov_gids_h = m_overlapped->get_dofs_gids().get_view<const gid_type*,Host>();

  // 1. Communicate GIDs lists to recv pids
  std::map<int,std::vector<gid_type>> send_pid2gids;
  for (int i=0; i<nimp; ++i) {
    auto pid = m_import_pids_h[i];
    auto lid = m_import_lids_h[i];
    send_pid2gids[pid].push_back(ov_gids_h[lid]);
  }
  for (auto& [pid,gids] : send_pid2gids) {
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend (gids.data(),gids.size(),mpi_gid_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::gather, creating send request (step 1)");
  }
  std::map<int,std::vector<gid_type>> recv_pid2gids;
  for (int i=0; i<nexp; ++i) {
    auto pid = m_export_pids_h[i];
    recv_pid2gids[pid].push_back(-1);
  }
  for (auto& [pid,gids] : recv_pid2gids) {
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv (gids.data(),gids.size(),mpi_gid_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::gather, creating recv request (step 1)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on send requests (step 1)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on recv requests (step 1)");
  send_req.clear();
  recv_req.clear();
  
  // 2. Communicate T's count for each GID to recv pids
  std::map<int,std::vector<int>> send_pid2count;
  for (int i=0; i<nimp; ++i) {
    auto pid = m_import_pids_h[i];
    auto lid = m_import_lids_h[i];
    send_pid2count[pid].push_back(src.at(lid).size());
  }
  for (auto& [pid,count] : send_pid2count) {
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend (count.data(),count.size(),MPI_INT,pid,tag,mpi_comm,&req),
                   "GridImportExport::gather, creating send request (step 2)");
  }

  std::map<int,std::vector<int>>      recv_pid2count;
  for (auto& [pid,gids] : recv_pid2gids) {
    recv_pid2count[pid].resize(gids.size());
  }
  for (auto& [pid,count] : recv_pid2count) {
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv (count.data(),count.size(),MPI_INT,pid,tag,mpi_comm,&req),
                   "GridImportExport::gather, creating recv request (step 2)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on send requests (step 2)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on recv requests (step 2)");
  send_req.clear();
  recv_req.clear();

  // 3. Pack and send the data.
  std::map<int,std::vector<T>> send_pid2data;
  const auto& ov_gid2lid = m_overlapped->get_gid2lid_map();
  for (const auto& [pid,gids] : send_pid2gids) {
    auto& data = send_pid2data[pid];
    for (auto g : gids) {
      auto lid = ov_gid2lid.at(g);
      data.insert(data.end(),src.at(lid).begin(),src.at(lid).end());
    }
  }
  for (auto& [pid,data] : send_pid2data) {
    auto& req = send_req.emplace_back();
    check_mpi_call(MPI_Isend (data.data(),data.size(),mpi_data_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::gather, creating send request (step 3)");
  }

  std::map<int,std::vector<T>> recv_pid2data;
  for (const auto& [pid,count] : recv_pid2count) {
    int n = 0;
    for (auto c : count) n += c;

    auto& data = recv_pid2data[pid];
    data.resize(n);
  }
  for (auto& [pid,data] : recv_pid2data) {
    auto& req = recv_req.emplace_back();
    check_mpi_call(MPI_Irecv (data.data(),data.size(),mpi_data_t,pid,tag,mpi_comm,&req),
                   "GridImportExport::gather, creating recv request (step 3)");
  }
  check_mpi_call(MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on send requests (step 3)");
  check_mpi_call(MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE),
                 "GridImportExport::gather, waiting on recv requests (step 3)");
  send_req.clear();
  recv_req.clear();

  // 4. Unpack received data in dst map
  const auto& recv_gid2lid = m_unique->get_gid2lid_map();
  for (const auto& [pid,data] : recv_pid2data) {
    const auto& count = recv_pid2count[pid];
    const auto& gids  = recv_pid2gids[pid];
    const int num_gids = count.size();
    for (int i=0, pos=0; i<num_gids; ++i) {
      auto lid = recv_gid2lid.at(gids[i]);
      auto curr_sz = dst[lid].size();
      dst[lid].resize(curr_sz+count[i]);
      for (int k=0; k<count[i]; ++k, ++pos) {
        dst[lid][curr_sz+k] = data[pos];
      }
    }
  }
}

} // namespace scream

#endif // EAMXX_GRID_IMPORT_EXPORT_DATA_HPP
