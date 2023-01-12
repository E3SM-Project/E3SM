/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Connectivity.hpp"
#include "ErrorDefs.hpp"

#include <array>
#include <algorithm>

namespace Homme
{

Connectivity::Connectivity ()
 : m_finalized    (false)
 , m_initialized  (false)
 , m_num_local_elements (-1)
 , m_max_corner_elements(-1)
{
  // Nothing to be done here
}

void Connectivity::set_comm (const Comm& comm)
{
  // Input comm must be valid (not storing a null MPI comm)
  assert (comm.mpi_comm()!=MPI_COMM_NULL);

  m_comm = comm;
}

void Connectivity::set_num_elements (const int num_local_elements)
{
  // We don't allow to change the number of elements once set. There may be downstream classes
  // that already read num_elements from this class, and would not be informed of the change.
  assert (!m_initialized);

  // Safety check
  assert (num_local_elements>=0);

  m_num_local_elements  = num_local_elements;

  m_num_connections = ExecViewManaged<int[NUM_CONNECTION_SHARINGS+1][NUM_CONNECTION_KINDS+1]>("Connections counts");
  h_num_connections = Kokkos::create_mirror_view(m_num_connections);

  // Initialize all counters to 0 (even the missing ones)
  Kokkos::deep_copy (h_num_connections,0);

  // Connectivity is now initialized
  m_initialized = true;
}

void Connectivity::set_max_corner_elements (const int max_corner_elements)
{
  m_max_corner_elements = max_corner_elements;
}

void Connectivity::add_connection (
  const int e1_lid, const int e1_gid, const std::uint8_t e1_dir, const std::uint8_t e1_dir_idx, const int e1_pid,
  const int e2_lid, const int e2_gid, const std::uint8_t e2_dir, const std::uint8_t e2_dir_idx, const int e2_pid)
{
  // Connectivity must be in initialized state but not in finalized state
  assert (m_initialized);
  assert (!m_finalized);

  // Comm must not be a null comm, otherwise checks on ranks may be misleading
  assert (m_comm.mpi_comm()!=MPI_COMM_NULL);

  // I believe edges appear twice in fortran, once per each ordering.
  // Here, we only need to store them once
  if (e1_pid == m_comm.rank()) {
    // There is no edge-to-corner connection. Either the elements share a corner or an edge.
    assert (m_helpers.CONNECTION_KIND[e1_dir]==m_helpers.CONNECTION_KIND[e2_dir]);

    // The elem lid on the local side must be valid!
    assert (e1_lid>=0);

    ConnectionInfo info;
    LidGidPos& local  = info.local;
    LidGidPos& remote = info.remote;

    // Local and remote elements lid and connection position
    local.lid  = e1_lid;
    local.gid  = e1_gid;
    local.dir  = e1_dir;
    local.dir_idx = e1_dir_idx;
    remote.lid = e2_lid;
    remote.gid = e2_gid;
    remote.dir = e2_dir;
    remote.dir_idx = e2_dir_idx;

    // Kind
    info.kind = etoi(m_helpers.CONNECTION_KIND[e1_dir]);

    // Direction
    static const Direction CONNECTION_DIRECTION[4][4] = {
      {Direction::BACKWARD, Direction::FORWARD , Direction::FORWARD,  Direction::BACKWARD},
      {Direction::FORWARD,  Direction::BACKWARD, Direction::BACKWARD, Direction::FORWARD },
      {Direction::FORWARD,  Direction::BACKWARD, Direction::BACKWARD, Direction::FORWARD },
      {Direction::BACKWARD, Direction::FORWARD , Direction::FORWARD,  Direction::BACKWARD}
    };
    info.direction = (local.dir < 4 ?
                      // edge
                      etoi(CONNECTION_DIRECTION[local.dir][remote.dir]) :
                      // corner
                      etoi(Direction::FORWARD));

    if (e2_pid != m_comm.rank()) {
      info.sharing = etoi(ConnectionSharing::SHARED);
      info.remote_pid = e2_pid;
    } else {
      info.remote_pid = -1;
      info.sharing = etoi(ConnectionSharing::LOCAL);
    }

    ++h_num_connections(info.sharing,info.kind);

    ucon_info.push_back(
      UConInfo{local.lid, local.gid, remote.lid, remote.gid,
               info.remote_pid,
               local.dir, local.dir_idx, remote.dir, remote.dir_idx,
               info.kind, info.sharing, info.direction});
  }
}

void Connectivity::finalize(const bool sanity_check)
{
  // Update counters for groups with same sharing/kind
  for (int kind=0; kind<NUM_CONNECTION_KINDS; ++kind) {
    for (int sharing=0; sharing<NUM_CONNECTION_SHARINGS; ++sharing) {
      h_num_connections(etoi(ConnectionSharing::ANY),kind) += h_num_connections(sharing,kind);
      h_num_connections(sharing,etoi(ConnectionKind::ANY)) += h_num_connections(sharing,kind);
    }
    h_num_connections(etoi(ConnectionKind::ANY),etoi(ConnectionKind::ANY)) += h_num_connections(etoi(ConnectionSharing::ANY),kind);
  }

  setup_ucon();

  m_finalized = true;
}

bool Connectivity::UConInfo::operator< (const UConInfo& o) const {
  // Sort on local (L/G)ID so that element data are contiguous.
  if (l_lid < o.l_lid) return true;
  if (l_lid > o.l_lid) return false;
  // Sort next on the direction so that corners having multiple connections have
  // these connections contiguous.
  if (l_dir < o.l_dir) return true;
  if (l_dir > o.l_dir) return false;
  // Sort the direction index so that we can run through them linearly with BFB
  // ordering.
  if (l_dir_idx < o.l_dir_idx) return true;
  if (l_dir_idx > o.l_dir_idx) return false;
  // These next are optional but might help a tiny bit with contiguity of memory
  // access.
  if (r_pid < o.r_pid) return true;
  if (r_pid > o.r_pid) return false;
  return r_gid < o.r_gid;
}

void Connectivity::setup_ucon () {
  const size_t nconn = ucon_info.size();

  d_ucon = decltype(d_ucon)("Unstructured Connections", nconn);
  const auto m_ucon = Kokkos::create_mirror_view(d_ucon);
  d_ucon_ptr = decltype(d_ucon_ptr)("Unstructured Connections Ptr",
                                    m_num_local_elements+1);
  h_ucon = decltype(h_ucon)("Unstructured Connections", nconn);
  h_ucon_ptr = Kokkos::create_mirror_view(d_ucon_ptr);

  h_ucon_ptr(0) = 0;
  h_ucon_ptr(m_num_local_elements) = nconn;
  if (nconn == 0) return;

  std::sort(ucon_info.begin(), ucon_info.end());

  { // Set up element pointers into ucon.
    int ie = 0;
    auto l_lid_curr = ucon_info[ie].l_lid;
    assert(l_lid_curr == ie);
    for (size_t i = 0; i < nconn; ++i) {
      assert(ie < m_num_local_elements);
      const auto& uci = ucon_info[i];
      if (uci.l_lid != l_lid_curr) {
        assert(l_lid_curr < uci.l_lid);
        l_lid_curr = uci.l_lid;
        ++ie;
        assert(uci.l_lid == ie);
        h_ucon_ptr(ie) = i;
      }
    }
    assert(ie == m_num_local_elements-1);
  }

  // Fill Info structs. h_ucon contains the large ConnectionInfo struct for use
  // in model initialization.
  for (size_t i = 0; i < nconn; ++i) {
    const auto& uci = ucon_info[i];
    auto& info = h_ucon(i);
    info.local.lid = uci.l_lid; info.local.gid = uci.l_gid;
    info.local.dir = uci.l_dir; info.local.dir_idx = uci.l_dir_idx;
    info.remote.lid = uci.r_lid; info.remote.gid = uci.r_gid;
    info.remote.dir = uci.r_dir; info.remote.dir_idx = uci.r_dir_idx;
    info.remote_pid = uci.r_pid;
    info.kind = uci.kind;
    info.sharing = uci.sharing;
    info.direction = uci.direction;
  }
  // m/d_ucon contain the much smaller HaloExchangeUnstructuredConnectionInfo
  // for use during halo exchanges.
  for (size_t i = 0; i < nconn; ++i) {
    const auto& uci = ucon_info[i];
    auto& info = m_ucon(i);
    info.local_lid = uci.l_lid;
    info.local_dir = uci.l_dir;
    info.kind = uci.kind;
    info.sharing = uci.sharing;
    info.direction = uci.direction;
    info.sharing_local_remote_iconn = -1;
    if (info.sharing == etoi(ConnectionSharing::LOCAL)) {
      const int je = h_ucon(i).remote.lid;
      const int jbeg = h_ucon_ptr(je), jend = h_ucon_ptr(je+1);
      for (int j = jbeg; j < jend; ++j)
        if (h_ucon(j).local.dir == h_ucon(i).remote.dir &&
            h_ucon(j).local.dir_idx == h_ucon(i).remote.dir_idx) {
          info.sharing_local_remote_iconn = j;
          break;
        }
      assert(info.sharing_local_remote_iconn >= 0);
    }
  }

  Kokkos::deep_copy(d_ucon, m_ucon);
  Kokkos::deep_copy(d_ucon_ptr, h_ucon_ptr);

  // Clear memory.
  ucon_info = decltype(ucon_info)();

  {
    /* Check that the connection counts are OK.
         Let max_elements_attached_to_node be the number of elements attached to
       an element corner. For example, in a regular planar mesh, it's 4. for RRM
       it's <= 7, as established in src/share/dimensions_mod.F90.
         Let max_corner_elem be the max number of elements attached to a corner
       that are not accounted for by the two edges and the parent element. Thus,
         max_corner_elem
             = max_elements_attached_to_node - (2 edge elements) - (1 parent element)
             = max_elements_attached_to_node - 3.
         Let max_neigh_edges be the max number of connections to an element:
           max_neigh_edges =   4*max_elements_attached_to_node
                             - 4*(1 edge element/corner to avoid duplication)
                             - 4*(1 parent element)
             = 4*max_elements_attached_to_node - 8.
    */
    const int max_elements_attached_to_node = 7;
    std::string msg;
    bool ok = true;
    if (m_max_corner_elements > max_elements_attached_to_node - 3) {
      msg = std::string("m_max_corner_elements is ")
        + std::to_string(m_max_corner_elements) + " but should be <= 7.\n";
      ok = false;
    }
    int max_conn_per_elem = 0; // our name for 'max_neigh_edges'
    for (int ie = 0; ie < m_num_local_elements; ++ie)
      max_conn_per_elem = std::max(max_conn_per_elem,
                                   h_ucon_ptr(ie+1) - h_ucon_ptr(ie));
    if (max_conn_per_elem > 4*max_elements_attached_to_node - 8) {
      const auto msg = std::string("max_conn_per_elem is ")
        + std::to_string(max_conn_per_elem) + " but should be <= 4*"
        + std::to_string(max_elements_attached_to_node) + " - 8.\n";
      ok = false;
    }
    for (int ie = 0; ie < m_num_local_elements; ++ie) {
      std::uint8_t max_corner_elem = 0;
      for (int k = h_ucon_ptr(ie); k < h_ucon_ptr(ie+1); ++k)
        max_corner_elem = std::max(max_corner_elem,
                                   h_ucon(k).local.dir_idx);
      if (max_corner_elem > m_max_corner_elements) {
        const auto msg = std::string("m_max_corner_elements is ")
          + std::to_string(m_max_corner_elements)
          + "but LID ie has " + std::to_string(static_cast<int>(max_corner_elem))
          + "\n";
        ok = false;
        break;
      }
    }
    if ( ! ok)
      Errors::runtime_abort(
        std::string("Connectivity::setup_ucon: At least one connection"
                    " count is wrong:\n") + msg);
  }

  {
    // Check that each element has four edge connections ordered S, N, W, E
    // before any of the corners. This assumption would have to be relaxed to
    // support domains with nonperiodic boundaries.
    //   Code that uses this assumption is tagged "assume:conn-edges-snwe".
    ConnectionHelpers h;
    for (int ie = 0; ie < m_num_local_elements; ++ie) {
      const int k = h_ucon_ptr(ie);
      bool ok = true;
      for (int i = 0; i < 4; ++i)
        ok = (ok &&
              h_ucon(k+i).kind == etoi(ConnectionKind::EDGE) &&
              h_ucon(k+i).local.dir == h.UNPACK_EDGES_ORDER[i] &&
              // In fact, S,N,W,E is numbered 0:3.
              h_ucon(k+i).local.dir == i);
      if ( ! ok)
        Errors::runtime_abort("Connectivity::setup_ucon: Element's first four"
                              " connections are not edges in S, N, W, E order.");
    }
  }
}

void Connectivity::clean_up()
{
  // Cleaning the elements counter
  Kokkos::deep_copy(h_num_connections,0);

  d_ucon = decltype(d_ucon)("", 0);
  h_ucon = decltype(h_ucon)("", 0);
  d_ucon_ptr = decltype(d_ucon_ptr)("", 0);
  h_ucon_ptr = decltype(h_ucon_ptr)("", 0);

  m_initialized = false;
  m_finalized   = false;
}

} // namespace Homme
