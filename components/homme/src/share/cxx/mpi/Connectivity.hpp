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
struct LidGidPos
{
  // The local and global IDs of the element.
  int lid, gid;
  // The cardinal direction (ConnectionName) and the index within the
  // connections to the element. For edges, dir == dir_idx in 0:3. For corners,
  // dir and dir_idx are not necessarily 1-1.
  std::uint8_t dir, dir_idx;
};

// A simple struct, storing a connection info. In addition to LidGidPos (on both local and
// remote element), it stores also whether the ordering is the same on both the element
// (relevant only for edge-type connections), and the process id of the remote element,
// which is only used if  the remote element is on a different process.
// Note: we store kind, sharing and direction already converted to ints
struct ConnectionInfo
{
  // This is only needed if the neighboring element is owned by a different process
  int remote_pid; // Process id owning the other side of the connection

  LidGidPos local;
  LidGidPos remote;

  std::uint8_t kind;     // etoi(ConnectionKind::EDGE)=0, etoi(ConnectionKind::CORNER)=1,  etoi(ConnectionSharing::MISSING)=2
  std::uint8_t sharing;  // etoi(ConnectionSharing::LOCAL)=0, etoi(ConnectionSharing::SHARED)=1, etoi(ConnectionSharing::MISSING)=2


  // The following is needed only for W/E/S/N edges, in case the ordering of the NP points is different in the two elements
  std::uint8_t direction;  //0=forward, 1=backward
};

// Just the data from the above that are needed on device during halo
// exchanges.
struct HaloExchangeUnstructuredConnectionInfo
{
  int local_lid;
  std::uint8_t local_dir, kind, sharing, direction;
  // If sharing == LOCAL, then record the connection index within the remote
  // element to which (dir, dir_idx) maps. This is for the pack phase's special
  // treatement of this case.
  int sharing_local_remote_iconn;
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
  void set_max_corner_elements (const int max_corner_elements);

  // An element's position is determined by
  // * its dir: 0-3: S, N, W, E edges; 4-7: corners and
  // * its order within the dir, the dir_idx, which is > 0 only for some
  //   corners in RRM grids.
  void add_connection (const int e1_lid, const int e1_gid, const std::uint8_t e1_dir, const std::uint8_t e1_dir_idx, const int e1_pid,
                       const int e2_lid, const int e2_gid, const std::uint8_t e2_dir, const std::uint8_t e2_dir_idx, const int e2_pid);

  void finalize (const bool sanity_check = true);

  void clean_up ();
  //@}

  //@name Getters
  //@{

  // Unstructured connections to handle RRM case. Connections for an element
  // having local ID ie are
  //   ucon(ucon_ptr(ie)):ucon(ucon_ptr(ie+1)-1).
  // Device d_ucon(i) is 1-1 with host h_ucon(i) but has a different, ~2.5x
  // smaller, format.
  ExecViewUnmanaged<const HaloExchangeUnstructuredConnectionInfo*> get_d_ucon () const { return d_ucon; }
  ExecViewUnmanaged<const int*> get_d_ucon_ptr () const { return d_ucon_ptr; }
  HostViewUnmanaged<const ConnectionInfo*> get_h_ucon () const { return h_ucon; }
  HostViewUnmanaged<const int*> get_h_ucon_ptr () const { return h_ucon_ptr; }

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
  int get_max_corner_elements    () const { return m_max_corner_elements; }

  bool is_initialized () const { return m_initialized; }
  bool is_finalized   () const { return m_finalized;   }

  const Comm& get_comm () const { return m_comm; }
  //@}

private:
  // An invalid id
  static constexpr int INVALID_ID = -1;
  // An invalid direction or direction_idx
  static constexpr std::uint8_t INVALID_DIR = 0xFF;

  Comm    m_comm;

  bool    m_finalized;
  bool    m_initialized;

  int     m_num_local_elements, m_max_corner_elements;

  ConnectionHelpers m_helpers;

  // TODO: do we need the counters on the device? It appears we never use them...
  ExecViewManaged<int[NUM_CONNECTION_SHARINGS+1][NUM_CONNECTION_KINDS+1]>             m_num_connections;
  ExecViewManaged<int[NUM_CONNECTION_SHARINGS+1][NUM_CONNECTION_KINDS+1]>::HostMirror h_num_connections;

  ExecViewManaged<HaloExchangeUnstructuredConnectionInfo*> d_ucon;
  ExecViewManaged<ConnectionInfo*>::HostMirror h_ucon;
  ExecViewManaged<int*>             d_ucon_ptr;
  ExecViewManaged<int*>::HostMirror h_ucon_ptr;
  ExecViewManaged<int*>             d_ucon_dir_ptr;
  ExecViewManaged<int*>::HostMirror h_ucon_dir_ptr;
  // Helper used to accumulated connections during add_connection phase. Emptied
  // in finalize. l_ is local; r_ is remote.
  struct UConInfo {
    int l_lid, l_gid, r_lid, r_gid, r_pid;
    std::uint8_t l_dir, l_dir_idx, r_dir, r_dir_idx, kind, sharing, direction;
    bool operator< (const UConInfo& o) const;
  };
  std::vector<UConInfo> ucon_info;
  // In finalize call, construct the unstructured connectivity data using
  // ucon_info.
  void setup_ucon();
};

} // namespace Homme

#endif // HOMMEXX_CONNECTIVITY_HPP
