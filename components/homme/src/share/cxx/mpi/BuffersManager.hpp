/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_MPI_BUFFERS_MANAGER_HPP
#define HOMMEXX_MPI_BUFFERS_MANAGER_HPP

#include "Types.hpp"

#include <vector>
#include <map>
#include <memory>

namespace Homme
{

// Forward declarations
class Connectivity;
class BoundaryExchange;

/*
 * BuffersManager: a class to handle the buffers needed by BoundaryExchange
 *
 * This class (BM) is responsible mainly for allocating buffers to be
 * used by the BoundaryExchange (BE) class. A single BM object can
 * be the buffers provider for several different BE objects, which
 * are referred to as 'customers'; in this case, the stored buffers
 * will have a size big enough to accommodate the most demanding needs
 * among those of all the customers. For more details about why one BM
 * may have more than one customer, see the explanation in the header
 * BoundaryExchange.hpp.
 *
 * The BM class stores 4 types of buffers, in the form of Kokkos View's:
 *
 *  - a local buffer: this can be non empty only if there is more than
 *    one element per processor. This buffer stores values that the BE
 *    customer needs to exchange with other elements on the same rank
 *  - a send and a recv buffer: these buffers are non empty only if
 *    there is more than one rank. These buffers will store the outgoing
 *    and incoming values that the BE customer is exchanging with
 *    other processes.
 *  - a blackhole send and recv buffer: these buffers are always used
 *    (regardless of the number of ranks), and they will be used to
 *    dump discarded values or read dummy zero values. Their need
 *    arises from the fact that there are missing connections between
 *    elements (see ConnectivityHelpers.hpp). In particular, there
 *    are some elements corners that are not shared with any other
 *    element per se. They are shared with other elements AS PART
 *    OF AN EDGE, but there is no other element that thouches the
 *    current one ONLY through that point. In order to avoid if's
 *    inside BE, we treat these corners as 'normal' corners, but
 *    the outgoing values will be stored in this blackhole send
 *    buffer (which is never read from), and the incoming values
 *    will be read from the blackhole recv buffer (which is filled
 *    with zeros).
 *  - an mpi send and recv buffer: these buffers are strictly linked
 *    to the send and recv ones. They indeed store the same pointers
 *    if the MPI can be performed using pointers in the Kokkos'
 *    Execution Space. This is always true for CPU/KNL builds, but
 *    it may or may not be true for GPU builds.
 *    The send/recv buffers are used to pack/unpack the data, while
 *    the mpi_send/mpi_recv buffers are used by MPI.
 *
 * The BM class also takes care of syncing the send/recv buffers
 * with the mpi_send/mpi_recv buffers, via a call to Kokkos::deep_copy,
 * which is a no-op if the MPIMemSpace=ExecMemSpace, that is, if
 * the MPI is performed using pointers on the Execution Space.
 *
 */

class BuffersManager
{
public:

  BuffersManager ();
  BuffersManager (std::shared_ptr<Connectivity> connectivity);
  ~BuffersManager ();

  // I'm not sure copying this class is a good idea.
  BuffersManager(const BuffersManager&) = delete;
  BuffersManager& operator= (const BuffersManager&) = delete;

  // Checks whether the connectivity is already set
  bool is_connectivity_set () const { return m_connectivity!=nullptr; }

  // Set the connectivity class
  void set_connectivity (std::shared_ptr<Connectivity> connectivity);

  // Ask the manager to re-check whether there is enough storage for all the BE's
  void check_for_reallocation ();

  // Check that the allocated views can handle the requested number of 2d/3d fields
  bool check_views_capacity (const int num_1d_fields, const int num_2d_fields,  const int num_3d_fields) const;

  // Allocate the buffers (overwriting possibly already allocated ones if needed)
  void allocate_buffers ();

  // Lock/unlock the buffers are busy
  void lock_buffers ();
  void unlock_buffers ();

  bool are_buffers_busy () const { return m_buffers_busy; }

  ExecViewUnmanaged<Real*> get_send_buffer           () const;
  ExecViewUnmanaged<Real*> get_recv_buffer           () const;
  ExecViewUnmanaged<Real*> get_local_buffer          () const;
  MPIViewUnmanaged<Real*>  get_mpi_send_buffer       () const;
  MPIViewUnmanaged<Real*>  get_mpi_recv_buffer       () const;
  ExecViewUnmanaged<Real*> get_blackhole_send_buffer () const;
  ExecViewUnmanaged<Real*> get_blackhole_recv_buffer () const;

  std::shared_ptr<Connectivity> get_connectivity () const { return m_connectivity; }

private:

  // Make BoundaryExchange a friend, so it can call the next four methods underneath
  friend class BoundaryExchange;

  // Adds/removes the given BoundaryExchange to/from the list of 'customers' of this class
  // Note: the only class that should call these methods is BoundaryExchange, so
  //       it can register/unregister itself as a customer
  void add_customer (BoundaryExchange* add_me);
  void remove_customer (BoundaryExchange* remove_me);
  // Deep copy the send/recv buffer to/from the mpi_send/recv buffer
  // Note: these are no-ops if MPIMemSpace=ExecMemSpace
  void sync_send_buffer (BoundaryExchange* customer);
  void sync_recv_buffer (BoundaryExchange* customer);

  // Small struct, to hold customer's needs. We could use an std::pair, but this is more verbose
  struct CustomerNeeds {
    size_t local_buffer_size;
    size_t mpi_buffer_size;

    bool operator== (const CustomerNeeds& rhs) {
      return local_buffer_size==rhs.local_buffer_size && mpi_buffer_size==rhs.mpi_buffer_size;
    }
  };

  // If necessary, updates buffers sizes so that there is enough storage to hold the required number of fields.
  // Note: this method does not (re)allocate views
  void update_requested_sizes (std::map<BoundaryExchange*,CustomerNeeds>::value_type& customer);

  // Computes the required storages
  void required_buffer_sizes (const int num_1d_fields, const int num_2d_fields, const int num_3d_fields,
                              size_t& mpi_buffer_size, size_t& local_buffer_size) const;

  // The number of customers
  size_t m_num_customers;

  // The sizes of the buffer
  size_t m_mpi_buffer_size;
  size_t m_local_buffer_size;

  // Used to check whether buffers are busy
  bool m_buffers_busy;

  // Used to check whether user can still request different sizes
  bool m_views_are_valid;

  // Customers of this BuffersManager, each with its local and mpi sizes
  std::map<BoundaryExchange*,CustomerNeeds>  m_customers;

  // The connectivity (needed to allocate buffers)
  std::shared_ptr<Connectivity> m_connectivity;

  // The buffers
  ExecViewManaged<Real*>  m_send_buffer;
  ExecViewManaged<Real*>  m_recv_buffer;
  ExecViewManaged<Real*>  m_local_buffer;

  // The mpi buffers (same as the previous send/recv buffers if MPIMemSpace=ExecMemSpace)
  MPIViewManaged<Real*>   m_mpi_send_buffer;
  MPIViewManaged<Real*>   m_mpi_recv_buffer;

  // The blackhole send/recv buffers (used for missing connections)
  ExecViewManaged<Real*>  m_blackhole_send_buffer;
  ExecViewManaged<Real*>  m_blackhole_recv_buffer;
};

inline void BuffersManager::sync_send_buffer (BoundaryExchange* customer)
{
  // Only customers can call this
  assert (m_customers.find(customer)!=m_customers.end());

  const size_t customer_mpi_buffer_size = m_customers.find(customer)->second.mpi_buffer_size;
  if (customer_mpi_buffer_size<m_mpi_buffer_size) {
    // Avoid copying more than we need
    MPIViewUnmanaged<Real*>  mpi_send_view(m_mpi_send_buffer.data(),customer_mpi_buffer_size);
    ExecViewUnmanaged<const Real*> send_view(m_send_buffer.data(),customer_mpi_buffer_size);
    Kokkos::deep_copy(mpi_send_view, send_view);
  } else {
    Kokkos::deep_copy(m_mpi_send_buffer, m_send_buffer);
  }
}

inline void BuffersManager::sync_recv_buffer (BoundaryExchange* customer)
{
  // Only customers can call this
  assert (m_customers.find(customer)!=m_customers.end());

  const size_t customer_mpi_buffer_size = m_customers.find(customer)->second.mpi_buffer_size;
  if (customer_mpi_buffer_size<m_mpi_buffer_size) {
    // Avoid copying more than we need
    MPIViewUnmanaged<const Real*>  mpi_recv_view(m_mpi_recv_buffer.data(),customer_mpi_buffer_size);
    ExecViewUnmanaged<Real*> recv_view(m_recv_buffer.data(),customer_mpi_buffer_size);
    Kokkos::deep_copy(recv_view, mpi_recv_view);
  } else {
    Kokkos::deep_copy(m_recv_buffer, m_mpi_recv_buffer);
  }
}

inline ExecViewUnmanaged<Real*>
BuffersManager::get_send_buffer () const
{
  // We ensure that the buffers are valid
  assert(m_views_are_valid);
  return m_send_buffer;
}

inline ExecViewUnmanaged<Real*>
BuffersManager::get_recv_buffer () const
{
  // We ensure that the buffers are valid
  assert(m_views_are_valid);
  return m_recv_buffer;
}

inline ExecViewUnmanaged<Real*>
BuffersManager::get_local_buffer () const
{
  // We ensure that the buffers are valid
  assert(m_views_are_valid);
  return m_local_buffer;
}

inline MPIViewUnmanaged<Real*>
BuffersManager::get_mpi_send_buffer() const
{
  // We ensure that the buffers are valid
  assert(m_views_are_valid);
  return m_mpi_send_buffer;
}

inline MPIViewUnmanaged<Real*>
BuffersManager::get_mpi_recv_buffer() const
{
  // We ensure that the buffers are valid
  assert(m_views_are_valid);
  return m_mpi_recv_buffer;
}

inline ExecViewUnmanaged<Real*>
BuffersManager::get_blackhole_send_buffer () const
{
  // We ensure that the buffers are valid
  assert(m_views_are_valid);
  return m_blackhole_send_buffer;
}

inline ExecViewUnmanaged<Real*>
BuffersManager::get_blackhole_recv_buffer () const
{
  // We ensure that the buffers are valid
  assert(m_views_are_valid);
  return m_blackhole_recv_buffer;
}

} // namespace Homme

#endif // HOMMEXX_MPI_BUFFERS_MANAGER_HPP
