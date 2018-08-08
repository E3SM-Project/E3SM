/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "BuffersManager.hpp"

#include "BoundaryExchange.hpp"
#include "Connectivity.hpp"

namespace Homme
{

BuffersManager::BuffersManager ()
 : m_num_customers     (0)
 , m_mpi_buffer_size   (0)
 , m_local_buffer_size (0)
 , m_buffers_busy      (false)
 , m_views_are_valid   (false)
{
  // The "fake" buffers used for MISSING connections. These do not depend on the requirements
  // from the custormers, so we can create them right away.
  constexpr size_t blackhole_buffer_size = 2 * NUM_LEV * VECTOR_SIZE;
  m_blackhole_send_buffer = ExecViewManaged<Real*>("blackhole array",blackhole_buffer_size);
  m_blackhole_recv_buffer = ExecViewManaged<Real*>("blackhole array",blackhole_buffer_size);
  Kokkos::deep_copy(m_blackhole_send_buffer,0.0);
  Kokkos::deep_copy(m_blackhole_recv_buffer,0.0);
}

BuffersManager::BuffersManager (std::shared_ptr<Connectivity> connectivity)
 : BuffersManager()
{
  set_connectivity(connectivity);
}

BuffersManager::~BuffersManager ()
{
  // Check that all the customers un-registered themselves
  assert (m_num_customers==0);

  // Check our buffers are not busy
  assert (!m_buffers_busy);
}

void BuffersManager::check_for_reallocation ()
{
  for (auto& it : m_customers) {
    update_requested_sizes (it);
  }
}

void BuffersManager::set_connectivity (std::shared_ptr<Connectivity> connectivity)
{
  // We don't allow a null connectivity, or a change of connectivity
  assert (connectivity && !m_connectivity);

  m_connectivity = connectivity;
}

bool BuffersManager::check_views_capacity (const int num_1d_fields, const int num_2d_fields, const int num_3d_fields) const
{
  size_t mpi_buffer_size, local_buffer_size;
  required_buffer_sizes (num_1d_fields, num_2d_fields, num_3d_fields, mpi_buffer_size, local_buffer_size);

  return (mpi_buffer_size<=m_mpi_buffer_size) &&
         (local_buffer_size<=m_local_buffer_size);
}

void BuffersManager::allocate_buffers ()
{
  // If views are marked as valid, they are already allocated, and no other
  // customer has requested a larger size
  if (m_views_are_valid) {
    return;
  }

  // The buffers used for packing/unpacking
  m_send_buffer  = ExecViewManaged<Real*>("send buffer",  m_mpi_buffer_size);
  m_recv_buffer  = ExecViewManaged<Real*>("recv buffer",  m_mpi_buffer_size);
  m_local_buffer = ExecViewManaged<Real*>("local buffer", m_local_buffer_size);

  // The buffers used in MPI calls
  m_mpi_send_buffer = Kokkos::create_mirror_view(decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  m_views_are_valid = true;

  // Tell to all our customers that they need to redo the setup of the internal buffer views
  for (auto& be_ptr : m_customers) {
    // Invalidate buffer views and requests in the customer (if none built yet, it's a no-op)
    be_ptr.first->clear_buffer_views_and_requests ();
  }
}

void BuffersManager::lock_buffers ()
{
  // Make sure we are not trying to lock buffers already locked
  assert (!m_buffers_busy);

  m_buffers_busy = true;
}

void BuffersManager::unlock_buffers ()
{
  // TODO: I am not checking if the buffers are locked. This allows to call
  //       the method twice in a row safely. Is this a bad idea?
  m_buffers_busy = false;
}

void BuffersManager::add_customer (BoundaryExchange* add_me)
{
  // We don't allow null customers (although this should never happen)
  assert (add_me!=nullptr);

  // We also don't allow re-registration
  assert (m_customers.find(add_me)==m_customers.end());

  // Add to the list of customers
  auto pair_it_bool = m_customers.emplace(add_me,CustomerNeeds{0,0});

  // Update the number of customers
  ++m_num_customers;

  // If this customer has already started the registration, we can already update the buffers sizes
  if (add_me->is_registration_started()) {
    update_requested_sizes(*pair_it_bool.first);
  }
}

void BuffersManager::remove_customer (BoundaryExchange* remove_me)
{
  // We don't allow null customers (although this should never happen)
  assert (remove_me!=nullptr);

  // Perhaps overzealous, but won't hurt: we should have customers
  assert (m_num_customers>0);

  // Find the customer
  auto it = m_customers.find(remove_me);

  // We don't allow removal of non-customers
  assert (it!=m_customers.end());

  // Remove the customer and its needs
  m_customers.erase(it);

  // Decrease number of customers
  --m_num_customers;
}

void BuffersManager::update_requested_sizes (typename std::map<BoundaryExchange*,CustomerNeeds>::value_type& customer)
{
  // Make sure connectivity is valid
  assert (m_connectivity && m_connectivity->is_finalized());

  // Make sure this is a customer
  assert (m_customers.find(customer.first)!=m_customers.end());

  // Get the number of fields that this customer has
  const int num_1d_fields = customer.first->get_num_1d_fields();
  const int num_2d_fields = customer.first->get_num_2d_fields();
  const int num_3d_fields = customer.first->get_num_3d_fields();

  // Compute the requested buffers sizes and compare with stored ones
  required_buffer_sizes (num_1d_fields, num_2d_fields, num_3d_fields, customer.second.mpi_buffer_size, customer.second.local_buffer_size);
  if (customer.second.mpi_buffer_size>m_mpi_buffer_size) {
    // Update the total
    m_mpi_buffer_size = customer.second.mpi_buffer_size;

    // Mark the views as invalid
    m_views_are_valid = false;
  }

  if(customer.second.local_buffer_size>m_local_buffer_size) {
    // Update the total
    m_local_buffer_size = customer.second.local_buffer_size;

    // Mark the views as invalid
    m_views_are_valid = false;
  }
}

void BuffersManager::required_buffer_sizes (const int num_1d_fields, const int num_2d_fields, const int num_3d_fields,
                                            size_t& mpi_buffer_size, size_t& local_buffer_size) const
{
  mpi_buffer_size = local_buffer_size = 0;

  // The buffer size for each connection kind
  // Note: for 2d/3d fields, we have 1 Real per GP (per level, in 3d). For 1d fields,
  //       we have 2 Real per level (max and min over element).
  int elem_buf_size[2];
  elem_buf_size[etoi(ConnectionKind::CORNER)] = num_1d_fields*2*NUM_LEV*VECTOR_SIZE + (num_2d_fields + num_3d_fields*NUM_LEV*VECTOR_SIZE) * 1;
  elem_buf_size[etoi(ConnectionKind::EDGE)]   = num_1d_fields*2*NUM_LEV*VECTOR_SIZE + (num_2d_fields + num_3d_fields*NUM_LEV*VECTOR_SIZE) * NP;

  // Compute the requested buffers sizes and compare with stored ones
  mpi_buffer_size += elem_buf_size[etoi(ConnectionKind::CORNER)] * m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::SHARED,ConnectionKind::CORNER);
  mpi_buffer_size += elem_buf_size[etoi(ConnectionKind::EDGE)]   * m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::SHARED,ConnectionKind::EDGE);

  local_buffer_size += elem_buf_size[etoi(ConnectionKind::CORNER)] * m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::LOCAL,ConnectionKind::CORNER);
  local_buffer_size += elem_buf_size[etoi(ConnectionKind::EDGE)]   * m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::LOCAL,ConnectionKind::EDGE);
}

} // namespace Homme
