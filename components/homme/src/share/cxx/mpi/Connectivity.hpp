/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CONNECTIVITY_HPP
#define HOMMEXX_CONNECTIVITY_HPP

#include "ConnectivityHelpers.hpp"

#include "Comm.hpp"

#include "Types.hpp"

namespace Homme
{
// A simple struct to store, for a connection between elements, the local/global id of the element
// and the position, meaning which neighbor this connection refers to (W/E/S/N/SW/SE/NW/NE)
// Much like std::tuple, but with more verbose members' names
struct LidGidPos
{
  int lid;
  int gid;
  int pos;
};

// An invalid id
constexpr int INVALID_ID = -1;

// A simple struct, storing a connection info. In addition to LidGidPos (on both local and
// remote element), it stores also whether the ordering is the same on both the element
// (relevant only for edge-type connections), and the process id of the remote element,
// which is only used if  the remote element is on a different process.
// Note: we store kind, sharing and direction already converted to ints
struct ConnectionInfo
{
  LidGidPos local;
  LidGidPos remote;

  int kind;     // etoi(ConnectionKind::EDGE)=0, etoi(ConnectionKind::CORNER)=1,  etoi(ConnectionSharing::MISSING)=2
  int sharing;  // etoi(ConnectionSharing::LOCAL)=0, etoi(ConnectionSharing::SHARED)=1, etoi(ConnectionSharing::MISSING)=2


  // The following is needed only for W/E/S/N edges, in case the ordering of the NP points is different in the two elements
  int direction;  //0=forward, 1=backward

  // This is only needed if the neighboring element is owned by a different process
  int remote_pid; // Process id owning the other side of the connection
};

// The connectivity class. It stores two lists of ConnectionInfo objects, one for
// local connections (both elements on process) and one for shared connections
// (one element is on a remote process). The latter require MPI work, while the
// former can be handled locally.
class Connectivity
{
public:

  Connectivity ();
  Connectivity& operator= (const Connectivity& src) = default;

  //@name Methods
  //@{

  void set_comm (const Comm& comm);

  void set_num_elements (const int num_local_elements);

  void add_connection (const int first_elem_lid,  const int first_elem_gid,  const int first_elem_pos,  const int first_elem_pid,
                       const int second_elem_lid, const int second_elem_gid, const int second_elem_pos, const int second_elem_pid);

  void finalize ();

  void clean_up ();
  //@}

  //@name Getters
  //@{

  // Get the view with all connections
  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<std::is_same<MemSpace,HostMemSpace>::value,
                 HostViewUnmanaged<const ConnectionInfo*[NUM_CONNECTIONS]>>::type
  get_connections () const { return h_connections; }

  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<std::is_same<MemSpace,ExecMemSpace>::value && !std::is_same<ExecMemSpace,HostMemSpace>::value,
                 ExecViewUnmanaged<const ConnectionInfo*[NUM_CONNECTIONS]>>::type
  get_connections () const { return m_connections; }

  // Get a particular connection
  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  const ConnectionInfo& get_connection (const int ie, const int iconn) const { return get_connections<MemSpace>()(ie,iconn); }

  // Get number of connections with given kind and sharing
  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<std::is_same<MemSpace,HostMemSpace>::value,int>::type
  get_num_connections (const ConnectionSharing sharing, const ConnectionKind kind) const { return h_num_connections(etoi(sharing), etoi(kind)); }

  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<std::is_same<MemSpace,ExecMemSpace>::value && !std::is_same<ExecMemSpace,HostMemSpace>::value,int>::type
  get_num_connections (const ConnectionSharing sharing, const ConnectionKind kind) const { return m_num_connections(etoi(sharing), etoi(kind)); }

  // Shortcuts of the previous getter for common sharing/kind pairs
  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  int get_num_connections        () const { return get_num_connections<MemSpace>(ConnectionSharing::ANY,   ConnectionKind::ANY); }
  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  int get_num_shared_connections () const { return get_num_connections<MemSpace>(ConnectionSharing::SHARED,ConnectionKind::ANY); }
  template<typename MemSpace>
  KOKKOS_INLINE_FUNCTION
  int get_num_local_connections  () const { return get_num_connections<MemSpace>(ConnectionSharing::LOCAL, ConnectionKind::ANY); }

  int get_num_local_elements     () const { return m_num_local_elements;  }

  bool is_initialized () const { return m_initialized; }
  bool is_finalized   () const { return m_finalized;   }

  const Comm& get_comm () const { return m_comm; }
  //@}

private:

  Comm    m_comm;

  bool    m_finalized;
  bool    m_initialized;

  int     m_num_local_elements;

  ConnectionHelpers m_helpers;

  // TODO: do we need the counters on the device? It appears we never use them...
  ExecViewManaged<int[NUM_CONNECTION_SHARINGS+1][NUM_CONNECTION_KINDS+1]>             m_num_connections;
  ExecViewManaged<int[NUM_CONNECTION_SHARINGS+1][NUM_CONNECTION_KINDS+1]>::HostMirror h_num_connections;

  ExecViewManaged<ConnectionInfo*[NUM_CONNECTIONS]>             m_connections;
  ExecViewManaged<ConnectionInfo*[NUM_CONNECTIONS]>::HostMirror h_connections;
};

} // namespace Homme

#endif // HOMMEXX_CONNECTIVITY_HPP
