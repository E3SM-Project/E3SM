/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Connectivity.hpp"

#include <array>
#include <algorithm>

namespace Homme
{

Connectivity::Connectivity ()
 : m_finalized    (false)
 , m_initialized  (false)
 , m_num_local_elements (-1)
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
  // Besides, does it really make sense? Does it add an interesting feature? I don't think so.
  assert (!m_initialized);

  // Safety check
  assert (num_local_elements>=0);

  m_num_local_elements  = num_local_elements;

  m_connections = ExecViewManaged<ConnectionInfo*[NUM_CONNECTIONS]>("Connections", m_num_local_elements);
  h_connections = Kokkos::create_mirror_view(m_connections);

  m_num_connections = ExecViewManaged<int[NUM_CONNECTION_SHARINGS+1][NUM_CONNECTION_KINDS+1]>("Connections counts");
  h_num_connections = Kokkos::create_mirror_view(m_num_connections);

  // Initialize all connections to MISSING
  // Note: we still include local element/position, since we need that even for
  //       missing connections! Everything else is set to invalid numbers (for safety)
  for (int ie=0; ie<m_num_local_elements; ++ie) {
    for (int iconn=0; iconn<NUM_CONNECTIONS; ++iconn) {
      ConnectionInfo& info = h_connections(ie,iconn);

      info.sharing = etoi(ConnectionSharing::MISSING);
      info.kind    = etoi(ConnectionKind::MISSING);

      info.local.lid = ie;
      info.local.pos = iconn;

      info.local.gid  = INVALID_ID;
      info.remote.lid = INVALID_ID;
      info.remote.gid = INVALID_ID;
      info.remote.pos = INVALID_ID;

    }
  }

  // Initialize all counters to 0 (even the missing ones)
  Kokkos::deep_copy (h_num_connections,0);

  // Connectivity is now initialized
  m_initialized = true;
}

void Connectivity::add_connection (const int first_elem_lid,  const int first_elem_gid,  const int first_elem_pos,  const int first_elem_pid,
                                   const int second_elem_lid, const int second_elem_gid, const int second_elem_pos, const int second_elem_pid)
{
  // Connectivity must be in initialized state but not in finalized state
  assert (m_initialized);
  assert (!m_finalized);

  // Comm must not be a null comm, otherwise checks on ranks may be misleading
  assert (m_comm.mpi_comm()!=MPI_COMM_NULL);

  // I believe edges appear twice in fortran, once per each ordering.
  // Here, we only need to store them once
  if (first_elem_pid==m_comm.rank())
  {
#ifndef NDEBUG
    // There is no edge-to-corner connection. Either the elements share a corner or an edge.
    assert (m_helpers.CONNECTION_KIND[first_elem_pos]==m_helpers.CONNECTION_KIND[second_elem_pos]);

    // The elem lid on the local side must be valid!
    assert (first_elem_lid>=0);
#endif

    ConnectionInfo& info = h_connections(first_elem_lid,first_elem_pos);
    LidGidPos& local  = info.local;
    LidGidPos& remote = info.remote;

    // Local and remote elements lid and connection position
    local.lid  = first_elem_lid;
    local.gid  = first_elem_gid;
    local.pos  = first_elem_pos;
    remote.lid = second_elem_lid;
    remote.gid = second_elem_gid;
    remote.pos = second_elem_pos;

    // Kind
    info.kind = etoi(m_helpers.CONNECTION_KIND[first_elem_pos]);

    // Direction
    info.direction = etoi(m_helpers.CONNECTION_DIRECTION[local.pos][remote.pos]);

    if (second_elem_pid!=m_comm.rank())
    {
      info.sharing = etoi(ConnectionSharing::SHARED);
      info.remote_pid = second_elem_pid;

    } else {
      info.remote_pid = -1;
      info.sharing = etoi(ConnectionSharing::LOCAL);
    }

    ++h_num_connections(info.sharing,info.kind);
  }
}

void Connectivity::finalize()
{
  // Sanity check: Homme does not allow less than 2*2 elements on each of the cube's faces, so each element
  // should have at most ONE missing connection
  constexpr int corners[NUM_CORNERS] = { etoi(ConnectionName::SWEST), etoi(ConnectionName::SEAST), etoi(ConnectionName::NWEST), etoi(ConnectionName::NEAST)};

  for (int ie=0; ie<m_num_local_elements; ++ie) {
#ifndef NDEBUG
    std::array<bool,NUM_CORNERS> missing = {{false, false, false, false}};
#endif

    for (int ic : corners) {
      if (h_connections(ie,ic).kind == etoi(ConnectionKind::MISSING)) {
#ifndef NDEBUG
        missing[ic % NUM_CORNERS] = true;
#endif

        // Just for tracking purposes
        //++h_num_connections(etoi(ConnectionSharing::MISSING),etoi(ConnectionKind::MISSING));
        h_num_connections(2,2) += 1;//etoi(ConnectionSharing::MISSING),etoi(ConnectionKind::MISSING)) += 1;
      }
    }
    assert (std::count(missing.cbegin(),missing.cend(),true)<=1);
  }

  // Updating counters for groups with same sharing/kind
  for (int kind=0; kind<NUM_CONNECTION_KINDS; ++kind) {
    for (int sharing=0; sharing<NUM_CONNECTION_SHARINGS; ++sharing) {
      h_num_connections(etoi(ConnectionSharing::ANY),kind) += h_num_connections(sharing,kind);
      h_num_connections(sharing,etoi(ConnectionKind::ANY)) += h_num_connections(sharing,kind);
    }
    h_num_connections(etoi(ConnectionKind::ANY),etoi(ConnectionKind::ANY)) += h_num_connections(etoi(ConnectionSharing::ANY),kind);
  }

  // Copying to device
  Kokkos::deep_copy(m_connections, h_connections);
  Kokkos::deep_copy(m_num_connections, h_num_connections);

  m_finalized = true;
}

void Connectivity::clean_up()
{
  m_connections = ExecViewManaged<ConnectionInfo*[NUM_CONNECTIONS]>("",0);
  Kokkos::deep_copy(m_num_connections,0);

  // Cleaning up also the host mirrors
  Kokkos::deep_copy(h_connections, m_connections);
  Kokkos::deep_copy(h_num_connections,0);

  // Cleaning the elements counter

  m_initialized = false;
  m_finalized   = false;
}

} // namespace Homme
