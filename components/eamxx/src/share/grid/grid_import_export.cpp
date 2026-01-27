#include "grid_import_export.hpp"

#include "share/field/field_utils.hpp"

namespace scream
{

GridImportExport::
GridImportExport (const std::shared_ptr<const AbstractGrid>& unique,
                  const std::shared_ptr<const AbstractGrid>& overlapped)
{
  EKAT_REQUIRE_MSG (unique!=nullptr, "Error! Input unique grid pointer is null.\n");
  EKAT_REQUIRE_MSG (overlapped!=nullptr, "Error! Input overlapped grid pointer is null.\n");

  EKAT_REQUIRE_MSG (unique->is_unique(),
      "Error! GridImportExport unique grid is not unique.\n");

  m_unique = unique;
  m_overlapped = overlapped;
  m_comm = unique->get_comm();

  // Note: we can't use the gids views from the grids, since we need to pass pointers
  //       to MPI bcast routines, which require pointers to nonconst data.
  //       Hence, create a clone, and grab a non-const pointer from it.
  const auto unique_gids_f = unique->get_dofs_gids().clone();
  const auto overlap_gids_f = overlapped->get_dofs_gids().clone();
  const auto    gids = unique_gids_f.get_view<gid_type*,Host>();
  const auto ov_gids = overlap_gids_f.get_view<gid_type*,Host>();

  int num_ov_gids = ov_gids.size();

  gid_type* data;
  std::vector<gid_type> pid_gids;
  std::map<int,std::vector<int>> pid2lids;

  // ------------------ Create import structures ----------------------- //

  // Resize output
  m_import_lids = decltype(m_import_lids)("",num_ov_gids);
  m_import_pids = decltype(m_import_pids)("",num_ov_gids);

  m_import_lids_h = Kokkos::create_mirror_view(m_import_lids);
  m_import_pids_h = Kokkos::create_mirror_view(m_import_pids);
  Kokkos::deep_copy(m_import_pids_h,-1);

  // We may have repeated gids. In that case, we want to update
  // the pids/lids arrays at all indices corresponding to the same gid
  // std::map<gid_type,std::vector<int>> gid2idx;
  // for (int i=0; i<num_ov_gids; ++i) {
  //   gid2idx[ov_gids[i]].push_back(i);
  // }
  auto ov_gid2lid = overlapped->get_gid2lid_map();

  // Let each rank bcast its src gids, so that other procs can
  // check against their dst grid
  int num_imports = 0;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    // Bcast dst gids count for this pid
    int num_gids_pid = gids.size();
    m_comm.broadcast(&num_gids_pid,1,pid);

    // Bcast gids
    if (pid==m_comm.rank()) {
      data = gids.data();
    } else {
      pid_gids.resize(num_gids_pid);
      data = pid_gids.data();
    }
    m_comm.broadcast(data,num_gids_pid,pid);

    // Checks if any of this pid's gids in the dst gids list
    for (int i=0; i<num_gids_pid; ++i) {
      auto it = ov_gid2lid.find(data[i]);
      if (it!=ov_gid2lid.end()) {
        pid2lids[pid].push_back(it->second);
        ++num_imports;
      }
    }
  }
  EKAT_REQUIRE_MSG (num_ov_gids==num_imports,
      "Error! Could not locate the owner of one of the dst grid GIDs.\n"
      "  - rank: " + std::to_string(m_comm.rank()) + "\n"
      "  - num found: " + std::to_string(num_imports) + "\n"
      "  - num dst gids: " + std::to_string(num_ov_gids) + "\n");
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    const auto& lids = pid2lids[pid];
    for (size_t i=0; i<lids.size(); ++i,++pos) {
      m_import_lids_h(pos) = lids[i];
      m_import_pids_h(pos) = pid;
    }
  }

  Kokkos::deep_copy(m_import_lids,m_import_lids_h);
  Kokkos::deep_copy(m_import_pids,m_import_pids_h);

  // ------------------ Create export structures ----------------------- //

  auto gid2lid = unique->get_gid2lid_map();

  // Let each rank bcast its src gids, so that other procs can
  // check against their dst grid
  // Note: we don't know a priori how many PIDs will need each
  // of our dofs, so we cannot insert in the export pids/lids views yet,
  // and must use a temporary map to store results
  int num_exports = 0;
  pid2lids.clear();
  for (int pid=0; pid<m_comm.size(); ++pid) {
    // Bcast dst gids count for this pid
    int num_gids_pid = ov_gids.size();
    m_comm.broadcast(&num_gids_pid,1,pid);

    // Bcast gids
    if (pid==m_comm.rank()) {
      data = ov_gids.data();
    } else {
      pid_gids.resize(num_gids_pid);
      data = pid_gids.data();
    }
    m_comm.broadcast(data,num_gids_pid,pid);

    // Checks if any of the src lids is in the list of needs of this pid
    for (int i=0; i<num_gids_pid; ++i) {
      auto it = gid2lid.find(data[i]);
      if (it!=gid2lid.end()) {
        pid2lids[pid].push_back(it->second);
        ++num_exports;
      }
    }

    // IMPORTANT! When building the import data, within each PID, we order
    // the list of imports according to the *remote* ordering. In order for
    // p2p messages to be consistent, the export data must order the
    // list of exports according to the *local* ordring.
    std::sort(pid2lids[pid].begin(),pid2lids[pid].end());
  }

  m_export_pids = view_1d<int>("",num_exports);
  m_export_lids = view_1d<int>("",num_exports);
  m_export_lids_h = Kokkos::create_mirror_view(m_export_lids);
  m_export_pids_h = Kokkos::create_mirror_view(m_export_pids);
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    const auto& lids = pid2lids[pid];
    for (size_t i=0; i<lids.size(); ++i,++pos) {
      m_export_lids_h(pos) = lids[i];
      m_export_pids_h(pos) = pid;
    }
  }

  // Kokkos::deep_copy(m_export_pids_count,pids_count_h);
  Kokkos::deep_copy(m_export_pids,m_export_pids_h);
  Kokkos::deep_copy(m_export_lids,m_export_lids_h);

  // Compute counts per pid
  m_num_exports_per_pid = view_1d<int>("",m_comm.size());
  m_num_exports_per_pid_h = Kokkos::create_mirror_view(m_num_exports_per_pid);
  for (size_t i=0; i<m_export_pids.size(); ++i) {
    ++m_num_exports_per_pid_h[m_export_pids_h[i]];
  }
  Kokkos::deep_copy(m_num_exports_per_pid,m_num_exports_per_pid_h);

  m_num_imports_per_pid = view_1d<int>("",m_comm.size());
  m_num_imports_per_pid_h = Kokkos::create_mirror_view(m_num_imports_per_pid);
  for (size_t i=0; i<m_import_pids.size(); ++i) {
    ++m_num_imports_per_pid_h[m_import_pids_h[i]];
  }
  Kokkos::deep_copy(m_num_imports_per_pid,m_num_imports_per_pid_h);
}

} // namespace scream
