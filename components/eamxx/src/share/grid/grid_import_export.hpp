#ifndef EAMXX_GRID_IMPORT_EXPORT_HPP
#define EAMXX_GRID_IMPORT_EXPORT_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/core/eamxx_types.hpp"       // For KokkosTypes

#include <ekat_comm.hpp>

#include <memory>

namespace scream
{

/*
 * Import/Export data is used to figure out where data
 * can be retrieved or sent (both in terms of remote
 * rank and local id) when transferring between
 * two grids. The terms import/export do not necessarily
 * imply that we are getting/sending data. Rather, they
 * have to do with which grid we're initiating the
 * transfer from. Namely,
 *   - import: the grid we have data on is potentially
 *             non unique, and the target grid is unique
 *   - export: the grid we have data on is unique, and
 *             the target grid is potentially non unique
 * The import/export data (pids/lids) can be used both on host
 * and device, for efficient pack/unpack methods.
 *
 * The data is organized as follows:
 *   - imports are indexed by the lids of the overlapped grid: import_lids(k) is the
 *     lid (in the overlapped grid) of the k-th import, and import_pids(k) is the rank
 *     owning the corresponding gid in the unique grid;
 *   - exports are indexed by the lids of the unique grid: export_lids(k) is the lid
 *     (in the unique grid) of the k-th export, and export_pids(k) is the rank that
 *     has that gid in its overlapped grid;
 *   - both are sorted by pid, and the offsets of each pid are stored too (view of size nranks+1);
 *   - within each pid, imports and exports are sorted by lid *in the unique grid*,
 *     so that if pid1 exports X to pid2, then the imports of pid2 from pid1 are ordered as
 *     the exports of pid1 to pid2. That is, p2p messages are consistent.
 * The construction is collective, and scalable (it uses a distributed directory of gids).
 */

class GridImportExport {
public:
  using KT = KokkosTypes<DefaultDevice>;
  using view_d = typename KT::view_1d<int>;
  using view_h = typename view_d::host_mirror_type;

  GridImportExport (const std::shared_ptr<const AbstractGrid>& unique,
                    const std::shared_ptr<const AbstractGrid>& overlapped);
  ~GridImportExport () = default;

  int num_exports () const { return m_num_exports; }
  int num_imports () const { return m_num_imports; }

  view_d export_pid_offsets () const { return m_export_pid_offset; }
  view_d import_pid_offsets () const { return m_import_pid_offset; }
  view_h export_pid_offsets_h () const { return m_export_pid_offset_h; }
  view_h import_pid_offsets_h () const { return m_import_pid_offset_h; }

  view_d import_pids () const { return m_import_pids; }
  view_d import_lids () const { return m_import_lids; }
  view_d export_pids () const { return m_export_pids; }
  view_d export_lids () const { return m_export_lids; }

  view_h export_lids_h () const { return m_export_lids_h; }

protected:

  using gid_type = AbstractGrid::gid_type;

  int m_num_exports = 0;
  int m_num_imports = 0;

  // The pids/lids of each import/export
  // All these arrays are sorted by pid. That is, all imports
  // for pid 1 come before imports for pid 2, and same for exports.
  view_d  m_import_pids;
  view_d  m_import_lids;
  view_d  m_export_pids;
  view_d  m_export_lids;
  view_h  m_export_lids_h;   // Only needed on host

  // Count how many import/exports per pid we have. These views have size nranks
  view_d  m_import_pid_offset;
  view_d  m_export_pid_offset;
  view_h  m_import_pid_offset_h;
  view_h  m_export_pid_offset_h;

  ekat::Comm    m_comm;
};

} // namespace scream

#endif // EAMXX_GRID_IMPORT_EXPORT_DATA_HPP
