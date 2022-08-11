/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "BoundaryExchange.hpp"

#include "MpiBuffersManager.hpp"
#include "KernelVariables.hpp"
#include "profiling.hpp"

#include "utilities/VectorUtils.hpp"

#define tstart(x)
#define tstop(x)

namespace Homme
{

void BoundaryExchange::build_buffer_views_and_requests()
{
  // If we already set the buffers before, then nothing to be done here
  if (m_buffer_views_and_requests_built) {
    return;
  }

  // Check that the MpiBuffersManager is present and was setup with enough storage
  assert (m_buffers_manager);

  // Ask the buffer manager to check for reallocation and then proceed with the allocation (if needed)
  // Note: if BM already knows about our needs, and buffers were already allocated, then
  //       these two calls should not change the internal state of the BM
  m_buffers_manager->check_for_reallocation();
  m_buffers_manager->allocate_buffers();

  // Note: this may look cryptic, so I'll try to explain what's about to happen.
  //       We want to set the send/recv buffers to point to:
  //         - a portion of send/recv_buffer if info.sharing=SHARED
  //         - a portion of local_buffer if info.sharing=LOCAL
  //         - the blackhole_send/recv if info.sharing=MISSING
  //       After reserving the buffer portion, update the offset by a given increment, depending on info.kind:
  //         - increment[CORNER]  = m_elem_buf_size[CORNER)] = 1  * (m_num_2d_fields + NUM_LEV*VECTOR_SIZE m_num_3d_fields)
  //         - increment[EDGE]    = m_elem_buf_size[EDGE)]   = NP * (m_num_2d_fields + NUM_LEV*VECTOR_SIZE m_num_3d_fields)
  //         - increment[MISSING] = 0 (point to the same blackhole)
  // Note: m_blackhole_send will be written many times, but will never be read from.
  //       Kind of like streaming to /dev/null. blackhole_recv will be read from sometimes
  //       (24 times, to be precise, one for each of the 3 corner connections on each of the
  //       cube's vertices), but it's never written into, so will always contain zeros (set by the constructor).

  HostViewManaged<size_t[3]> h_buf_offset("");
  Kokkos::deep_copy(h_buf_offset, 0);

  // The amount of Real's used in a connection on a single level:
  //  - 2d/3d field exchange 1 Real per GP
  //  - 1d fields exchange 2 Real per level (max/min over element)
  HostViewManaged<int[3]> h_increment_1d("increment_1d");
  HostViewManaged<int[3]> h_increment_2d("increment_2d");
  h_increment_1d[etoi(ConnectionKind::EDGE)]      =  2;
  h_increment_1d[etoi(ConnectionKind::CORNER)]    =  2;
  h_increment_1d[etoi(ConnectionKind::MISSING)]   =  0;
  h_increment_2d[etoi(ConnectionKind::EDGE)]    = NP;
  h_increment_2d[etoi(ConnectionKind::CORNER)]  =  1;
  h_increment_2d[etoi(ConnectionKind::MISSING)] =  0;
  HostViewManaged<int[3]> h_increment_3d = h_increment_2d;

  // Since we access the manager many times, we may as well call lock once and store the shared_ptr.
  auto buffers_manager = m_buffers_manager;

  using local_buf_ptr_type = decltype( buffers_manager->get_local_buffer().data());
  using local_buf_val_type = decltype(*buffers_manager->get_local_buffer().data());
  HostViewManaged<Pointer<local_buf_ptr_type, local_buf_val_type>[3]> h_all_send_buffers("");
  HostViewManaged<Pointer<local_buf_ptr_type, local_buf_val_type>[3]> h_all_recv_buffers("");

  h_all_send_buffers[etoi(ConnectionSharing::LOCAL)]   = buffers_manager->get_local_buffer().data();
  h_all_send_buffers[etoi(ConnectionSharing::SHARED)]  = buffers_manager->get_send_buffer().data();
  h_all_send_buffers[etoi(ConnectionSharing::MISSING)] = buffers_manager->get_blackhole_send_buffer().data();
  h_all_recv_buffers[etoi(ConnectionSharing::LOCAL)]   = buffers_manager->get_local_buffer().data();
  h_all_recv_buffers[etoi(ConnectionSharing::SHARED)]  = buffers_manager->get_recv_buffer().data();
  h_all_recv_buffers[etoi(ConnectionSharing::MISSING)] = buffers_manager->get_blackhole_recv_buffer().data();

  // Create buffer views
  m_send_1d_buffers = decltype(m_send_1d_buffers)("1d send buffer", m_num_elems, m_num_1d_fields);
  m_recv_1d_buffers = decltype(m_recv_1d_buffers)("1d recv buffer", m_num_elems, m_num_1d_fields);
  m_send_2d_buffers = decltype(m_send_2d_buffers)("2d send buffer", m_num_elems, m_num_2d_fields);
  m_recv_2d_buffers = decltype(m_recv_2d_buffers)("2d recv buffer", m_num_elems, m_num_2d_fields);
  m_send_3d_buffers = decltype(m_send_3d_buffers)("3d send buffer", m_num_elems, m_num_3d_fields);
  m_recv_3d_buffers = decltype(m_recv_3d_buffers)("3d recv buffer", m_num_elems, m_num_3d_fields);
  m_send_3d_int_buffers = decltype(m_send_3d_int_buffers)("3d interface send buffer", m_num_elems, m_num_3d_int_fields);
  m_recv_3d_int_buffers = decltype(m_recv_3d_int_buffers)("3d interface recv buffer", m_num_elems, m_num_3d_int_fields);

  std::vector<int> slot_idx_to_elem_conn_pair, pids, pid_offsets;
  init_slot_idx_to_elem_conn_pair(slot_idx_to_elem_conn_pair, pids, pid_offsets);

  // NOTE: I wanted to do this setup in parallel, on the execution space, but there
  //       is a reduction hidden. In particular, we need to access buf_offset atomically,
  //       so that it is not update while we are still using it. One solution would be to
  //       use Kokkos::atomic_fetch_add, but it may kill the concurrency. And given that
  //       we do not have a ton of concurrency in this setup phase, and given that it is
  //       precisely only a setup phase, we may as well do things serially on the host,
  //       then deep_copy back to device
  ConnectionHelpers helpers;
  auto h_send_1d_buffers = Kokkos::create_mirror_view(m_send_1d_buffers);
  auto h_recv_1d_buffers = Kokkos::create_mirror_view(m_recv_1d_buffers);
  auto h_send_2d_buffers = Kokkos::create_mirror_view(m_send_2d_buffers);
  auto h_recv_2d_buffers = Kokkos::create_mirror_view(m_recv_2d_buffers);
  auto h_send_3d_buffers = Kokkos::create_mirror_view(m_send_3d_buffers);
  auto h_recv_3d_buffers = Kokkos::create_mirror_view(m_recv_3d_buffers);
  auto h_send_3d_int_buffers = Kokkos::create_mirror_view(m_send_3d_int_buffers);
  auto h_recv_3d_int_buffers = Kokkos::create_mirror_view(m_recv_3d_int_buffers);
  auto h_connections = m_connectivity->get_connections<HostMemSpace>();
  for (int k = 0; k < m_num_elems*NUM_CONNECTIONS; ++k) {
    const int ie = slot_idx_to_elem_conn_pair[k] / NUM_CONNECTIONS;
    const int iconn = slot_idx_to_elem_conn_pair[k] % NUM_CONNECTIONS;
    {
      const ConnectionInfo& info = h_connections(ie, iconn);

      const LidGidPos local = info.local;

      auto send_buffer = h_all_send_buffers[info.sharing];
      auto recv_buffer = h_all_recv_buffers[info.sharing];

      for (int ifield=0; ifield<m_num_1d_fields; ++ifield) {
        h_send_1d_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Scalar[2][NUM_LEV]>(
          reinterpret_cast<Scalar*>(send_buffer.get() + h_buf_offset[info.sharing]));
        h_recv_1d_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Scalar[2][NUM_LEV]>(
          reinterpret_cast<Scalar*>(recv_buffer.get() + h_buf_offset[info.sharing]));
        h_buf_offset[info.sharing] += h_increment_1d[info.kind]*NUM_LEV*VECTOR_SIZE;
      }
      for (int ifield=0; ifield<m_num_2d_fields; ++ifield) {
        h_send_2d_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Real*>(
          send_buffer.get() + h_buf_offset[info.sharing], helpers.CONNECTION_SIZE[info.kind]);
        h_recv_2d_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Real*>(
          recv_buffer.get() + h_buf_offset[info.sharing], helpers.CONNECTION_SIZE[info.kind]);
        h_buf_offset[info.sharing] += h_increment_2d[info.kind];
      }
      for (int ifield=0; ifield<m_num_3d_fields; ++ifield) {
        h_send_3d_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Scalar*[NUM_LEV]>(
          reinterpret_cast<Scalar*>(send_buffer.get() + h_buf_offset[info.sharing]),
          helpers.CONNECTION_SIZE[info.kind]);
        h_recv_3d_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Scalar*[NUM_LEV]>(
          reinterpret_cast<Scalar*>(recv_buffer.get() + h_buf_offset[info.sharing]),
          helpers.CONNECTION_SIZE[info.kind]);
        h_buf_offset[info.sharing] += h_increment_3d[info.kind]*NUM_LEV*VECTOR_SIZE;
      }
      for (int ifield=0; ifield<m_num_3d_int_fields; ++ifield) {
        h_send_3d_int_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Scalar*[NUM_LEV_P]>(
          reinterpret_cast<Scalar*>(send_buffer.get() + h_buf_offset[info.sharing]),
          helpers.CONNECTION_SIZE[info.kind]);
        h_recv_3d_int_buffers(local.lid, ifield, local.pos) = ExecViewUnmanaged<Scalar*[NUM_LEV_P]>(
          reinterpret_cast<Scalar*>(recv_buffer.get() + h_buf_offset[info.sharing]),
          helpers.CONNECTION_SIZE[info.kind]);
        h_buf_offset[info.sharing] += h_increment_3d[info.kind]*NUM_LEV_P*VECTOR_SIZE;
      }
    }
  }
  Kokkos::deep_copy(m_send_1d_buffers, h_send_1d_buffers);
  Kokkos::deep_copy(m_recv_1d_buffers, h_recv_1d_buffers);
  Kokkos::deep_copy(m_send_2d_buffers, h_send_2d_buffers);
  Kokkos::deep_copy(m_recv_2d_buffers, h_recv_2d_buffers);
  Kokkos::deep_copy(m_send_3d_buffers, h_send_3d_buffers);
  Kokkos::deep_copy(m_recv_3d_buffers, h_recv_3d_buffers);
  Kokkos::deep_copy(m_send_3d_int_buffers, h_send_3d_int_buffers);
  Kokkos::deep_copy(m_recv_3d_int_buffers, h_recv_3d_int_buffers);

#ifndef NDEBUG
  // Sanity check: compute the buffers sizes for this boundary exchange, and checking that the final offsets match them
  size_t mpi_buffer_size = 0;
  size_t local_buffer_size = 0;

  mpi_buffer_size   += m_elem_buf_size[etoi(ConnectionKind::CORNER)] *
    m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::SHARED, ConnectionKind::CORNER);
  mpi_buffer_size   += m_elem_buf_size[etoi(ConnectionKind::EDGE)]   *
    m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::SHARED, ConnectionKind::EDGE);

  local_buffer_size += m_elem_buf_size[etoi(ConnectionKind::CORNER)] *
    m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::LOCAL, ConnectionKind::CORNER);
  local_buffer_size += m_elem_buf_size[etoi(ConnectionKind::EDGE)]   *
    m_connectivity->get_num_connections<HostMemSpace>(ConnectionSharing::LOCAL, ConnectionKind::EDGE);

  assert (h_buf_offset[etoi(ConnectionSharing::LOCAL)]==local_buffer_size);
  assert (h_buf_offset[etoi(ConnectionSharing::SHARED)]==mpi_buffer_size);
#endif // NDEBUG

  {
    auto mpi_comm = m_connectivity->get_comm().mpi_comm();
    const auto& connections = m_connectivity->get_connections<HostMemSpace>();
    const int npids = pids.size();
    free_requests();
    m_send_requests.resize(npids);
    m_recv_requests.resize(npids);
    MPIViewManaged<Real*>::pointer_type send_ptr = buffers_manager->get_mpi_send_buffer().data();
    MPIViewManaged<Real*>::pointer_type recv_ptr = buffers_manager->get_mpi_recv_buffer().data();
    int offset = 0;
    for (int ip = 0; ip < npids; ++ip) {
      int count = 0;
      for (int k = pid_offsets[ip]; k < pid_offsets[ip+1]; ++k) {
        const int ie = slot_idx_to_elem_conn_pair[k] / NUM_CONNECTIONS;
        const int iconn = slot_idx_to_elem_conn_pair[k] % NUM_CONNECTIONS;
        const ConnectionInfo& info = connections(ie, iconn);
        count += m_elem_buf_size[info.kind];
      }
      HOMMEXX_MPI_CHECK_ERROR(MPI_Send_init(send_ptr + offset, count, MPI_DOUBLE,
                                            pids[ip], m_exchange_type, mpi_comm,
                                            &m_send_requests[ip]),
                              m_connectivity->get_comm().mpi_comm());
      HOMMEXX_MPI_CHECK_ERROR(MPI_Recv_init(recv_ptr + offset, count, MPI_DOUBLE,
                                            pids[ip], m_exchange_type, mpi_comm,
                                            &m_recv_requests[ip]),
                              m_connectivity->get_comm().mpi_comm());
      offset += count;
    }
  }

  // Now the buffer views and the requests are built
  m_buffer_views_and_requests_built = true;
}

} // namespace Homme
