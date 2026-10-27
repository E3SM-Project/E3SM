#include "grid_import_export.hpp"

#include "share/util/gid_directory.hpp"
#include "share/util/eamxx_alltoallv.hpp"
#include "share/field/field_utils.hpp"

#include <algorithm>
#include <numeric>

namespace scream
{

GridImportExport::
GridImportExport (const std::shared_ptr<const AbstractGrid>& unique,
                  const std::shared_ptr<const AbstractGrid>& overlapped)
{
  EKAT_REQUIRE_MSG (unique!=nullptr, "Error! Input unique grid pointer is null.\n");
  EKAT_REQUIRE_MSG (overlapped!=nullptr, "Error! Input overlapped grid pointer is null.\n");

  m_comm = unique->get_comm();

  const int nranks = m_comm.size();

  auto unique_gids_h = unique->get_dofs_gids().get_view<const gid_type*,Host>();
  auto ov_gids_h     = overlapped->get_dofs_gids().get_view<const gid_type*,Host>();

  const int num_ov_gids = ov_gids_h.size();

  // The idea is to avoid every rank looking at the gids of all other ranks (which is not
  // scalable). Instead, we use a distributed directory, which tells us who owns a gid.
  //  1. Register the gids of the unique grid in the directory (this also checks that they are unique)
  //  2. Ask the directory who owns each gid in the overlapped grid (and what its lid is)
  //     This gives us the imports.
  //  3. Tell each owner which of its lids we want from it. This gives the owners their exports.

  // ------------------ Find the owners ----------------------- //

  GidDirectory<gid_type> directory (m_comm,unique_gids_h);
  EKAT_REQUIRE_MSG (not directory.has_duplicates(),
      "Error! GridImportExport unique grid is not unique.\n");

  const auto owners = directory.lookup(ov_gids_h);
  const auto& owner_pids = owners.pids;
  const auto& owner_lids = owners.lids;
  for (int i=0; i<num_ov_gids; ++i) {
    EKAT_REQUIRE_MSG (owner_pids[i]>=0,
        "Error! Could not locate the owner of one of the dst grid GIDs.\n"
        "  - rank: " + std::to_string(m_comm.rank()) + "\n"
        "  - gid : " + std::to_string(ov_gids_h[i]) + "\n"
        "  - num dst gids: " + std::to_string(num_ov_gids) + "\n");
  }

  // ------------------ Create import structures ----------------------- //

  // Sort the imports by pid. Within each pid, sort them by the *remote* lid. The owners
  // will sort their exports in the same way, so that p2p messages are consistent.
  m_num_imports = num_ov_gids;
  std::vector<int> order (num_ov_gids);
  std::iota(order.begin(),order.end(),0);
  std::sort(order.begin(),order.end(),
            [&](const int a, const int b) {
              return owner_pids[a]<owner_pids[b] or
                    (owner_pids[a]==owner_pids[b] and (owner_lids[a]<owner_lids[b] or
                    (owner_lids[a]==owner_lids[b] and a<b)));
            });

  m_import_lids = decltype(m_import_lids)("",m_num_imports);
  m_import_pids = decltype(m_import_pids)("",m_num_imports);
  auto import_lids_h = Kokkos::create_mirror_view(m_import_lids);
  auto import_pids_h = Kokkos::create_mirror_view(m_import_pids);

  std::vector<int> imp_count (nranks,0);
  std::vector<int> requested_lids (m_num_imports);
  for (int k=0; k<m_num_imports; ++k) {
    const int i = order[k];
    import_lids_h(k) = i;
    import_pids_h(k) = owner_pids[i];
    requested_lids[k]  = owner_lids[i];
    ++imp_count[owner_pids[i]];
  }
  Kokkos::deep_copy(m_import_lids,import_lids_h);
  Kokkos::deep_copy(m_import_pids,import_pids_h);

  // ------------------ Create export structures ----------------------- //

  // Send to each owner the lids we want from it. What we receive are the lids that others want from us,
  // already sorted by requesting pid, and (within each pid) by lid.
  std::vector<int> exported_lids, exp_count;
  alltoallv(m_comm,requested_lids,imp_count,exported_lids,exp_count,false);

  m_num_exports = exported_lids.size();
  m_export_pids = view_d("",m_num_exports);
  m_export_lids = view_d("",m_num_exports);
  m_export_lids_h = Kokkos::create_mirror_view(m_export_lids);
  auto export_pids_h = Kokkos::create_mirror_view(m_export_pids);
  for (int pid=0, pos=0; pid<nranks; ++pid) {
    for (int k=0; k<exp_count[pid]; ++k, ++pos) {
      m_export_lids_h(pos) = exported_lids[pos];
      export_pids_h(pos) = pid;
    }
  }
  Kokkos::deep_copy(m_export_pids,export_pids_h);
  Kokkos::deep_copy(m_export_lids,m_export_lids_h);

  // ------------------ Create offset views ----------------------- //

  // Compute offsets of each pid in import/export views
  m_export_pid_offset = view_d("",nranks+1);
  m_export_pid_offset_h = Kokkos::create_mirror_view(m_export_pid_offset);
  m_export_pid_offset_h[0] = 0;
  for (int i=0; i<nranks; ++i) {
    m_export_pid_offset_h[i+1] = m_export_pid_offset_h[i]+exp_count[i];
  }
  Kokkos::deep_copy(m_export_pid_offset,m_export_pid_offset_h);

  m_import_pid_offset = view_d("",nranks+1);
  m_import_pid_offset_h = Kokkos::create_mirror_view(m_import_pid_offset);
  m_import_pid_offset_h[0] = 0;
  for (int i=0; i<nranks; ++i) {
    m_import_pid_offset_h[i+1] = m_import_pid_offset_h[i]+imp_count[i];
  }
  Kokkos::deep_copy(m_import_pid_offset,m_import_pid_offset_h);
}

} // namespace scream
