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
  int pos;

  // ------------------ Create import structures ----------------------- //

  // Resize output
  m_import_lids = decltype(m_import_lids)("",num_ov_gids);
  m_import_pids = decltype(m_import_pids)("",num_ov_gids);

  m_import_lids_h = Kokkos::create_mirror_view(m_import_lids);
  m_import_pids_h = Kokkos::create_mirror_view(m_import_pids);
  Kokkos::deep_copy(m_import_pids_h,-1);

  // We may have repeated gids. In that case, we want to update
  // the pids/lids arrays at all indices corresponding to the same gid
  auto ov_gid2lid = overlapped->get_gid2lid_map();

  // Let each rank bcast its src gids, so that other procs can
  // check against their dst grid
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
        ++m_num_imports;
      }
    }
  }
  EKAT_REQUIRE_MSG (num_ov_gids==m_num_imports,
      "Error! Could not locate the owner of one of the dst grid GIDs.\n"
      "  - rank: " + std::to_string(m_comm.rank()) + "\n"
      "  - num found: " + std::to_string(m_num_imports) + "\n"
      "  - num dst gids: " + std::to_string(num_ov_gids) + "\n");

  pos=0;
  for (const auto& [pid,lids] : pid2lids) {
    for (auto lid : lids) {
      m_import_lids_h(pos) = lid;
      m_import_pids_h(pos) = pid;
      ++pos;
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
        ++m_num_exports;
      }
    }

    // IMPORTANT! When building the import data, within each PID, we order
    // the list of imports according to the *remote* ordering. In order for
    // p2p messages to be consistent, the export data must order the
    // list of exports according to the *local* ordring.
    std::sort(pid2lids[pid].begin(),pid2lids[pid].end());
  }

  m_export_pids = view_d("",m_num_exports);
  m_export_lids = view_d("",m_num_exports);
  m_export_lids_h = Kokkos::create_mirror_view(m_export_lids);
  m_export_pids_h = Kokkos::create_mirror_view(m_export_pids);

  pos=0;
  for (const auto& [pid,lids] : pid2lids) {
    for (auto lid : lids) {
      m_export_lids_h(pos) = lid;
      m_export_pids_h(pos) = pid;
      ++pos;
    }
  }

  Kokkos::deep_copy(m_export_pids,m_export_pids_h);
  Kokkos::deep_copy(m_export_lids,m_export_lids_h);

  // ------------------ Create offset views ----------------------- //

  // Compute offsets of each pid in import/export views
  m_export_pid_offset = view_d("",m_comm.size()+1);
  m_export_pid_offset_h = Kokkos::create_mirror_view(m_export_pid_offset);
  std::vector<int> exp_count (m_comm.size(),0);
  for (int i=0; i<m_num_exports; ++i) {
    ++exp_count[m_export_pids_h[i]];
  }
  m_export_pid_offset_h[0] = 0;
  for (int i=0; i<m_comm.size(); ++i) {
    m_export_pid_offset_h[i+1] = m_export_pid_offset_h[i]+exp_count[i];
  }
  Kokkos::deep_copy(m_export_pid_offset,m_export_pid_offset_h);

  m_import_pid_offset = view_d("",m_comm.size()+1);
  m_import_pid_offset_h = Kokkos::create_mirror_view(m_import_pid_offset);
  std::vector<int> imp_count (m_comm.size(),0);
  for (int i=0; i<m_num_imports; ++i) {
    ++imp_count[m_import_pids_h[i]];
  }
  m_import_pid_offset_h[0] = 0;
  for (int i=0; i<m_comm.size(); ++i) {
    m_import_pid_offset_h[i+1] = m_import_pid_offset_h[i]+imp_count[i];
  }
  Kokkos::deep_copy(m_import_pid_offset,m_import_pid_offset_h);
}

} // namespace scream
