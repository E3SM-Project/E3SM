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

// ======================== IMPLEMENTATION ======================== //

// Separating these allocations into a small routine works around a Cuda 10/GCC
// 7/C++14 internal error.
template <typename A, typename B>
static void alloc3d (A& a, B& b, int ne, int a2, int b2) {
  a = ExecViewManaged<ExecViewManaged<Scalar[NP][NP][NUM_LEV]>**>("3d fields", ne, a2);
  b = ExecViewManaged<ExecViewManaged<Scalar[NP][NP][NUM_LEV_P]>**>("3d interface fields", ne, b2);
}

BoundaryExchange::BoundaryExchange()
{
  m_num_1d_fields = 0;
  m_num_2d_fields = 0;
  m_num_3d_fields = 0;
  m_num_3d_int_fields = 0;

  m_connectivity    = std::shared_ptr<Connectivity>();
  m_buffers_manager = std::shared_ptr<MpiBuffersManager>();

  m_num_elems = -1;

  // Prohibit registration until the number of fields has been set
  m_registration_started   = false;
  m_registration_completed = false;

  // There is no buffer view or request yet
  m_buffer_views_and_requests_built = false;

  // We start with a clean class
  m_cleaned_up = true;
  m_send_pending = false;
  m_recv_pending = false;
}

BoundaryExchange::BoundaryExchange(std::shared_ptr<Connectivity> connectivity, std::shared_ptr<MpiBuffersManager> buffers_manager)
 : BoundaryExchange ()
{
  // Set the connectivity
  set_connectivity (connectivity);

  // Set the buffers manager
  set_buffers_manager (buffers_manager);
}

BoundaryExchange::~BoundaryExchange()
{
  clean_up ();

  // It may be that we never really used this object, and never set the BM...
  if (m_buffers_manager) {
    // Remove me as a customer of the BM
    m_buffers_manager->remove_customer(this);
  }
}

void BoundaryExchange::set_connectivity (std::shared_ptr<Connectivity> connectivity)
{
  // Functionality only available before registration starts
  assert (!m_registration_started && !m_registration_completed);

  // Make sure it is a valid connectivity (does not need to be initialized/finalized yet)
  // Also, replacing the connectivity could have unintended side-effects; better prohibit it.
  // Besides, when can it be useful?
  assert (connectivity && !m_connectivity);

  // If the buffers manager is set and it stores a connectivity, it must match the input one
  assert (!m_buffers_manager || m_buffers_manager->get_connectivity()==connectivity);

  // Set the connectivity
  m_connectivity = connectivity;
  m_num_elems = connectivity->get_num_local_elements();
}

void BoundaryExchange::set_buffers_manager (std::shared_ptr<MpiBuffersManager> buffers_manager)
{
  // Functionality available only before the registration is completed
  assert (!m_registration_completed);

  // Make sure it is a valid pointer. Also, replacing the buffers manager
  // could have unintended side-effects; better prohibit it.
  // Besides, when can it be useful?
  assert (buffers_manager && !m_buffers_manager);

  // If the buffers manager stores a connectivity, and we already have one set, they must match
  assert (!buffers_manager->is_connectivity_set() || !(m_connectivity) || buffers_manager->get_connectivity()==m_connectivity);

  // Set the internal pointer
  m_buffers_manager = buffers_manager;

  // Set the connectivity in the buffers manager, if not already set
  if (!m_buffers_manager->is_connectivity_set() && m_connectivity) {
    m_buffers_manager->set_connectivity(m_connectivity);
  }

  // If I don't store a connectivity, take it from the buffers manager (if it has one)
  if (m_buffers_manager->is_connectivity_set() && !m_connectivity) {
    set_connectivity(m_buffers_manager->get_connectivity());
  }

  // Add myself as a customer of the BM
  m_buffers_manager->add_customer(this);
}

void BoundaryExchange::set_num_fields (const int num_1d_fields, const int num_2d_fields, const int num_3d_fields, const int num_3d_int_fields)
{
  // We don't allow to call this method twice in a row. If you want to change the number of fields,
  // you need to call clean_up first, to get a fresh new BoundaryExchange.
  // Note: if you do call clean_up and then again set_num_field and the new number of fields
  //       are smaller than the previous ones, MpiBuffersManager will NOT shrink the buffers, so
  //       you may be left with buffers that are larger than what you need.
  assert (m_cleaned_up);

  // Make sure the connectivity is valid: must at least be a valid pointer, and be initialized, i.e.,
  // store a valid number of elements, but may be finalized later (before registration_completed call though)
  assert (m_connectivity && m_connectivity->is_initialized());

  // We strongly advocate for not using the same BE object for both 'standard' exchange and min/max exchange
  assert (!(num_1d_fields>0 && (num_2d_fields>0 || num_3d_fields>0 || num_3d_int_fields>0)));

  // Note: we do not set m_num_1d_fields, m_num_2d_fields and m_num_3d_fields, since we will use them as
  //       progressive indices while adding fields. Then, during registration_completed,
  //       we will check that they match the 2nd dimension of m_1d_fields, m_2d_fields and m_3d_fields.

  // Create the fields views
  m_1d_fields = decltype(m_1d_fields)("1d fields", m_num_elems, num_1d_fields);
  m_2d_fields = decltype(m_2d_fields)("2d fields", m_num_elems, num_2d_fields);
  if (NUM_LEV==NUM_LEV_P) {
    // If NUM_LEV=NUM_LEVP, we can use the same 3d buffers for both midpoints and interface quantities
    alloc3d(m_3d_fields, m_3d_int_fields, m_num_elems, num_3d_fields+num_3d_int_fields, 0);
  } else {
    alloc3d(m_3d_fields, m_3d_int_fields, m_num_elems, num_3d_fields, num_3d_int_fields);
  }

  // Now we can start register fields
  m_registration_started   = true;
  m_registration_completed = false;

  // We're not all clean
  m_cleaned_up = false;
}

void BoundaryExchange::clean_up()
{
  if (m_cleaned_up) {
    // Perhaps not possible, but just in case
    return;
  }

  // Check that we are not still transmitting
  assert (!m_send_pending && !m_recv_pending);

  // Clear stored fields
  m_1d_fields = decltype(m_1d_fields)("m_1d_fields", 0, 0);
  m_2d_fields = decltype(m_2d_fields)("m_2d_fields", 0, 0);
  alloc3d(m_3d_fields, m_3d_int_fields, 0, 0, 0);

  m_num_1d_fields = 0;
  m_num_2d_fields = 0;
  m_num_3d_fields = 0;
  m_num_3d_int_fields = 0;

  // If we clean up, we need to reset the number of fields
  m_registration_started   = false;
  m_registration_completed = false;

  // Clean buffer views and requests
  clear_buffer_views_and_requests();

  // Now we're all cleaned
  m_cleaned_up = true;
}

void BoundaryExchange::registration_completed()
{
  // If everything is already set up, just return
  if (m_registration_completed) {
    // TODO: should we prohibit two consecutive calls of this method? It seems harmless, so I'm allowing it
    return;
  }

  // TODO: should we assert that m_registration_started=true? Or simply return if not? Can calling this
  //       method without a call to registration started be dangerous? Not sure...

  // At this point, the connectivity MUST be finalized already, and the buffers manager must be set already
  assert (m_connectivity && m_connectivity->is_finalized());
  assert (m_buffers_manager);

  // Create the MPI data types, for corners and edges
  // Note: this is the size per element, per connection. It is the number of Real's to send/receive to/from the neighbor
  // Note: for 2d/3d fields, we have 1 Real per GP (per level, in 3d). For 1d fields,
  //       we have 2 Real per level (max and min over element).

  const int single_ptr_buf_size = m_num_2d_fields + m_num_3d_fields*NUM_LEV*VECTOR_SIZE + m_num_3d_int_fields*NUM_LEV_P*VECTOR_SIZE;
  m_elem_buf_size[etoi(ConnectionKind::CORNER)] = m_num_1d_fields*2*NUM_LEV*VECTOR_SIZE + single_ptr_buf_size * 1;
  m_elem_buf_size[etoi(ConnectionKind::EDGE)]   = m_num_1d_fields*2*NUM_LEV*VECTOR_SIZE + single_ptr_buf_size * NP;

  // Determine what kind of BE is this (exchange or exchange_min_max)
  m_exchange_type = m_num_1d_fields>0 ? MPI_EXCHANGE_MIN_MAX : MPI_EXCHANGE;

  // Prohibit further registration of fields, and allow exchange
  m_registration_started   = false;
  m_registration_completed = true;

  // Optimistically build buffers here. If registration is called with largest
  // BufferManager user first, then building will occur just once, in the
  // prim_init2 call.
  build_buffer_views_and_requests();
}

void BoundaryExchange::exchange () {
  exchange(nullptr);
}

void BoundaryExchange::exchange (ExecViewUnmanaged<const Real * [NP][NP]> rspheremp) {
  exchange(&rspheremp);
}

void BoundaryExchange::exchange (const ExecViewUnmanaged<const Real * [NP][NP]>* rspheremp)
{
  // Check that the registration has completed first
  assert (m_registration_completed);

  // Check that this object is setup to perform exchange and not exchange_min_max
  assert (m_exchange_type==MPI_EXCHANGE);

  // I am not sure why and if we could have this scenario, but just in case. I think MPI *may* go bananas in this case
  if (m_num_2d_fields+m_num_3d_fields+m_num_3d_int_fields==0) {
    return;
  }

  // If this is the first time we call the exchange method, or if the MpiBuffersManager has performed a reallocation
  // since the last time this method was called, we need to rebuild all our internal buffer views
  if (!m_buffer_views_and_requests_built) {
    build_buffer_views_and_requests();
  }

  // Hey, if some process can already send me stuff while I'm still packing, that's ok
  if ( ! m_recv_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Startall(m_recv_requests.size(), m_recv_requests.data()),
                            m_connectivity->get_comm().mpi_comm());
  m_recv_pending = true;

  // ---- Pack and send ---- //
  pack_and_send ();

  // --- Recv and unpack --- //
  recv_and_unpack (rspheremp);
}

void BoundaryExchange::exchange_min_max ()
{
  // Check that the registration has completed first
  assert (m_registration_completed);

  // Check that this object is setup to perform exchange_min_max and not exchange
  assert (m_exchange_type==MPI_EXCHANGE_MIN_MAX);

  // I am not sure why and if we could have this scenario, but just in case. I think MPI *may* go bananas in this case
  if (m_num_1d_fields==0) {
    return;
  }

  // If this is the first time we call the exchange method, or if the MpiBuffersManager has performed a reallocation
  // since the last time this method was called, we need to rebuild all our internal buffer views
  if (!m_buffer_views_and_requests_built) {
    build_buffer_views_and_requests();
  }

  // Hey, if some process can already send me stuff while I'm still packing, that's ok
  if ( ! m_recv_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Startall(m_recv_requests.size(), m_recv_requests.data()),
                            m_connectivity->get_comm().mpi_comm());
  m_recv_pending = true;

  // ---- Pack and send ---- //
  pack_and_send_min_max ();

  // --- Recv and unpack --- //
  recv_and_unpack_min_max ();
}

void BoundaryExchange::pack_and_send ()
{
  tstart("be pack_and_send");
  // The registration MUST be completed by now
  // Note: this also implies connectivity and buffers manager are valid
  assert (m_registration_completed);

  // Check that this object is setup to perform exchange and not exchange_min_max
  assert (m_exchange_type==MPI_EXCHANGE);

  // I am not sure why and if we could have this scenario, but just in case. I think MPI *may* go bananas in this case
  if (m_num_2d_fields+m_num_3d_fields+m_num_3d_int_fields==0) {
    return;
  }

  // Check that buffers are not locked by someone else, then lock them
  assert (!m_buffers_manager->are_buffers_busy());
  m_buffers_manager->lock_buffers();

  // If this is the first time we call this method, or if the MpiBuffersManager has performed a reallocation
  // since the last time this method was called, AND we are calling this method manually, without relying
  // on the exchange method to call it, then we need to rebuild all our internal buffer views
  if (!m_buffer_views_and_requests_built) {
    tstart("be build_buffer_views_and_requests");
    build_buffer_views_and_requests();
    tstop("be build_buffer_views_and_requests");
  }

  // ---- Pack ---- //
  // First, pack 2d fields (if any)...
  auto connections = m_connectivity->get_connections<ExecMemSpace>();
  if (m_num_2d_fields>0) {
    auto fields_2d = m_2d_fields;
    auto send_2d_buffers = m_send_2d_buffers;
    const ConnectionHelpers helpers;
    Kokkos::parallel_for(MDRangePolicy<ExecSpace, 3>({0, 0, 0}, {m_num_elems, NUM_CONNECTIONS, m_num_2d_fields}, {1, 1, 1}),
                         KOKKOS_LAMBDA(const int ie, const int iconn, const int ifield) {
      const ConnectionInfo& info = connections(ie, iconn);
      const LidGidPos& field_lidpos  = info.local;
      // For the buffer, in case of local connection, use remote info. In fact, while with shared connections the
      // mpi call will take care of "copying" data to the remote recv buffer in the correct remote element lid,
      // for local connections we need to manually copy on the remote element lid. We can do it here
      const LidGidPos& buffer_lidpos = info.sharing==etoi(ConnectionSharing::LOCAL) ? info.remote : info.local;

      // Note: if it is an edge and the remote edge is in the reverse order, we read the field_lidpos points backwards
      const auto& pts = helpers.CONNECTION_PTS[info.direction][field_lidpos.pos];
      for (int k=0; k<helpers.CONNECTION_SIZE[info.kind]; ++k) {
        send_2d_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos)(k) = fields_2d(field_lidpos.lid, ifield)(pts[k].ip, pts[k].jp);
      }
    });
  }
  // ...then pack 3d fields (if any)...
  if (m_num_3d_fields>0) {
    auto fields_3d = m_3d_fields;
    auto send_3d_buffers = m_send_3d_buffers;
    const auto num_3d_fields = m_num_3d_fields;
    if (OnGpu<ExecSpace>::value) {
      const ConnectionHelpers helpers;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*m_num_3d_fields*NUM_CONNECTIONS*NUM_LEV),
        KOKKOS_LAMBDA(const int it) {
          const int ie = it / (num_3d_fields*NUM_CONNECTIONS*NUM_LEV);
          const int ifield = (it / (NUM_CONNECTIONS*NUM_LEV)) % num_3d_fields;
          const int iconn = (it / NUM_LEV) % NUM_CONNECTIONS;
          const int ilev = it % NUM_LEV;
          const ConnectionInfo& info = connections(ie, iconn);
          const LidGidPos& field_lidpos = info.local;
          // For the buffer, in case of local connection, use remote info. In fact, while with shared connections the
          // mpi call will take care of "copying" data to the remote recv buffer in the correct remote element lid,
          // for local connections we need to manually copy on the remote element lid. We can do it here
          const LidGidPos& buffer_lidpos = info.sharing==etoi(ConnectionSharing::LOCAL) ? info.remote : info.local;

          // Note: if it is an edge and the remote edge is in the reverse order, we read the field_lidpos points backwards
          const auto& pts = helpers.CONNECTION_PTS[info.direction][field_lidpos.pos];
          const auto& sb = send_3d_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos);
          const auto& f3 = fields_3d(field_lidpos.lid, ifield);
          for (int k=0; k<helpers.CONNECTION_SIZE[info.kind]; ++k) {
            sb(k, ilev) = f3(pts[k].ip, pts[k].jp, ilev);
          }
        });
    } else {
      const auto num_parallel_iterations = m_num_elems*m_num_3d_fields;
      ThreadPreferences tp;
      tp.max_threads_usable = NUM_CONNECTIONS;
      tp.max_vectors_usable = NUM_LEV;
      const auto threads_vectors =
        DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
          num_parallel_iterations, tp);
      const auto policy = Kokkos::TeamPolicy<ExecSpace>(
        num_parallel_iterations, threads_vectors.first, threads_vectors.second);
      HOMMEXX_STATIC const ConnectionHelpers helpers;
      Kokkos::parallel_for(
        policy,
        KOKKOS_LAMBDA(const TeamMember& team) {
          Homme::KernelVariables kv(team, num_3d_fields);
          const int ie = kv.ie;
          const int ifield = kv.iq;
          for (int iconn = 0; iconn < 8; ++iconn) {
            const ConnectionInfo& info = connections(ie, iconn);
            if (info.kind == etoi(ConnectionSharing::MISSING)) continue;
            const LidGidPos& field_lidpos = info.local;
            const LidGidPos& buffer_lidpos = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                              info.remote :
                                              info.local);
            const auto& pts = helpers.CONNECTION_PTS[info.direction][field_lidpos.pos];
            const auto& sb = send_3d_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos);
            const auto& f3 = fields_3d(field_lidpos.lid, ifield);
            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(kv.team, helpers.CONNECTION_SIZE[info.kind]),
              [&] (const int& k) {
                const auto& ip = pts[k].ip;
                const auto& jp = pts[k].jp;
                auto* const sbp = &sb(k, 0);
                const auto* const f3p = &f3(ip, jp, 0);
                Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                  [&] (const int& ilev) {
                    sbp[ilev] = f3p[ilev];
                  });
              });
          }
        });
    }
  }
  // ...then pack 3d interface fields (if any)
  if (m_num_3d_int_fields>0) {
    auto fields = m_3d_int_fields;
    auto send_buffers = m_send_3d_int_buffers;
    const auto num_fields = m_num_3d_int_fields;
    if (OnGpu<ExecSpace>::value) {
      const ConnectionHelpers helpers;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*num_fields*NUM_CONNECTIONS*NUM_LEV_P),
        KOKKOS_LAMBDA(const int it) {
          const int ie = it / (num_fields*NUM_CONNECTIONS*NUM_LEV_P);
          const int ifield = (it / (NUM_CONNECTIONS*NUM_LEV_P)) % num_fields;
          const int iconn = (it / NUM_LEV_P) % NUM_CONNECTIONS;
          const int ilev = it % NUM_LEV_P;
          const ConnectionInfo& info = connections(ie, iconn);
          const LidGidPos& field_lidpos = info.local;
          // For the buffer, in case of local connection, use remote info. In fact, while with shared connections the
          // mpi call will take care of "copying" data to the remote recv buffer in the correct remote element lid,
          // for local connections we need to manually copy on the remote element lid. We can do it here
          const LidGidPos& buffer_lidpos = info.sharing==etoi(ConnectionSharing::LOCAL) ? info.remote : info.local;

          // Note: if it is an edge and the remote edge is in the reverse order, we read the field_lidpos points backwards
          const auto& pts = helpers.CONNECTION_PTS[info.direction][field_lidpos.pos];
          const auto& sb = send_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos);
          const auto& f = fields(field_lidpos.lid, ifield);
          for (int k=0; k<helpers.CONNECTION_SIZE[info.kind]; ++k) {
            sb(k, ilev) = f(pts[k].ip, pts[k].jp, ilev);
          }
        });
    } else {
      const auto num_parallel_iterations = m_num_elems*num_fields;
      ThreadPreferences tp;
      tp.max_threads_usable = NUM_CONNECTIONS;
      tp.max_vectors_usable = NUM_LEV_P;
      const auto threads_vectors =
        DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
          num_parallel_iterations, tp);
      const auto policy = Kokkos::TeamPolicy<ExecSpace>(
        num_parallel_iterations, threads_vectors.first, threads_vectors.second);
      HOMMEXX_STATIC const ConnectionHelpers helpers;
      Kokkos::parallel_for(
        policy,
        KOKKOS_LAMBDA(const TeamMember& team) {
          Homme::KernelVariables kv(team, num_fields);
          const int ie = kv.ie;
          const int ifield = kv.iq;
          for (int iconn = 0; iconn < 8; ++iconn) {
            const ConnectionInfo& info = connections(ie, iconn);
            if (info.kind == etoi(ConnectionSharing::MISSING)) continue;
            const LidGidPos& field_lidpos = info.local;
            const LidGidPos& buffer_lidpos = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                              info.remote :
                                              info.local);
            const auto& pts = helpers.CONNECTION_PTS[info.direction][field_lidpos.pos];
            const auto& sb = send_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos);
            const auto& f = fields(field_lidpos.lid, ifield);
            Kokkos::parallel_for(
              Kokkos::TeamThreadRange(kv.team, helpers.CONNECTION_SIZE[info.kind]),
              [&] (const int& k) {
                const auto& ip = pts[k].ip;
                const auto& jp = pts[k].jp;
                auto* const sbp = &sb(k, 0);
                const auto* const fp = &f(ip, jp, 0);
                Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(kv.team, NUM_LEV_P),
                  [&] (const int& ilev) {
                    sbp[ilev] = fp[ilev];
                  });
              });
          }
        });
    }
  }
  ExecSpace::impl_static_fence();

  // ---- Send ---- //
  tstart("be sync_send_buffer");
  m_buffers_manager->sync_send_buffer(this); // Deep copy send_buffer into mpi_send_buffer (no op if MPI is on device)
  tstop("be sync_send_buffer");
  tstart("be send");
  if ( ! m_send_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Startall(m_send_requests.size(), m_send_requests.data()),
                            m_connectivity->get_comm().mpi_comm());

  // Notify a send is ongoing
  m_send_pending = true;
  tstop("be pack_and_send");
}

void BoundaryExchange::recv_and_unpack () {
  recv_and_unpack(nullptr);
}

void BoundaryExchange::recv_and_unpack (const ExecViewUnmanaged<const Real * [NP][NP]>* rspheremp)
{
  tstart("be recv_and_unpack");
  tstart("be recv_and_unpack book");
  // The registration MUST be completed by now
  // Note: this also implies connectivity and buffers manager are valid
  assert (m_registration_completed);

  // Check that this object is setup to perform exchange and not exchange_min_max
  assert (m_exchange_type==MPI_EXCHANGE);

  // I am not sure why and if we could have this scenario, but just in case. I
  // think MPI *may* go bananas in this case
  if (m_num_2d_fields+m_num_3d_fields==0) {
    return;
  }

  // If I am doing pack_and_send and recv_and_unpack manually (rather than
  // through 'exchange'), then I need to start receiving now (otherwise it is
  // done already inside 'exchange')
  if (!m_recv_pending) {
    // If you are doing send/recv manually, don't call recv without a send, or
    // else you'll be stuck waiting later on
    assert (m_send_pending);

    if ( ! m_recv_requests.empty())
      HOMMEXX_MPI_CHECK_ERROR(MPI_Startall(m_recv_requests.size(), m_recv_requests.data()),
                              m_connectivity->get_comm().mpi_comm());
    m_recv_pending = true;
  }
  tstop("be recv_and_unpack book");

  // ---- Recv ---- //
  tstart("be recv waitall");
  if ( ! m_recv_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Waitall(m_recv_requests.size(), m_recv_requests.data(), MPI_STATUSES_IGNORE),
                            m_connectivity->get_comm().mpi_comm()); // Wait for all data to arrive
  m_recv_pending = false;
  tstop("be recv waitall");

  tstart("be recv_and_unpack book");
  m_buffers_manager->sync_recv_buffer(this);

  tstop("be recv_and_unpack book");

  // --- Unpack --- //
  // First, unpack 2d fields (if any)...
  if (m_num_2d_fields>0) {
    auto fields_2d = m_2d_fields;
    auto recv_2d_buffers = m_recv_2d_buffers;
    const ConnectionHelpers helpers;
    Kokkos::parallel_for(MDRangePolicy<ExecSpace, 2>({0, 0}, {m_num_elems, m_num_2d_fields}, {1, 1}),
                         KOKKOS_LAMBDA(const int ie, const int ifield) {
      for (int k=0; k<NP; ++k) {
        for (int iedge : helpers.UNPACK_EDGES_ORDER) {
          fields_2d(ie, ifield)(helpers.CONNECTION_PTS_FWD[iedge][k].ip,
                                helpers.CONNECTION_PTS_FWD[iedge][k].jp)
            += recv_2d_buffers(ie, ifield, iedge)[k];
        }
      }
      for (int icorner : helpers.UNPACK_CORNERS_ORDER) {
        if (recv_2d_buffers(ie, ifield, icorner).size() > 0) {
          fields_2d(ie, ifield)(helpers.CONNECTION_PTS_FWD[icorner][0].ip,
                                helpers.CONNECTION_PTS_FWD[icorner][0].jp)
            += recv_2d_buffers(ie, ifield, icorner)[0];
        }
      }
    });
  }
  // ...then unpack 3d fields (if any)...
  if (m_num_3d_fields>0) {
    auto fields_3d = m_3d_fields;
    auto recv_3d_buffers = m_recv_3d_buffers;
    const auto num_3d_fields = m_num_3d_fields;
    if (OnGpu<ExecSpace>::value) {
      const ConnectionHelpers helpers;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*m_num_3d_fields*NUM_LEV),
        KOKKOS_LAMBDA(const int it) {
          const int ie = it / (num_3d_fields*NUM_LEV);
          const int ifield = (it / NUM_LEV) % num_3d_fields;
          const int ilev = it % NUM_LEV;
          const auto& f3 = fields_3d(ie, ifield);
          for (int k=0; k<NP; ++k) {
            for (int iedge : helpers.UNPACK_EDGES_ORDER) {
              f3(helpers.CONNECTION_PTS_FWD[iedge][k].ip,
                 helpers.CONNECTION_PTS_FWD[iedge][k].jp, ilev)
                += recv_3d_buffers(ie, ifield, iedge)(k, ilev);
            }
          }
          for (int icorner : helpers.UNPACK_CORNERS_ORDER) {
            if (recv_3d_buffers(ie, ifield, icorner).size() > 0)
              f3(helpers.CONNECTION_PTS_FWD[icorner][0].ip,
                 helpers.CONNECTION_PTS_FWD[icorner][0].jp, ilev)
                += recv_3d_buffers(ie, ifield, icorner)(0, ilev);
          }
        });
      if (rspheremp) {
        const auto rsmp = *rspheremp;
        Kokkos::parallel_for(
          Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*m_num_3d_fields*NP*NP*NUM_LEV),
          KOKKOS_LAMBDA(const int it) {
            const int ie = it / (num_3d_fields*NUM_LEV*NP*NP);
            const int ifield = (it / (NP*NP*NUM_LEV)) % num_3d_fields;
            const int i = (it / (NP*NUM_LEV)) % NP;
            const int j = (it / NUM_LEV) % NP;
            const int ilev = it % NUM_LEV;
            fields_3d(ie, ifield)(i, j, ilev) *= rsmp(ie, i, j);
          });
      }
    } else {
      const auto num_parallel_iterations = m_num_elems*m_num_3d_fields;
      Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(num_parallel_iterations, 1, NUM_LEV),
        KOKKOS_LAMBDA(const TeamMember& team) {
          Homme::KernelVariables kv(team, num_3d_fields);
          const int ie = kv.ie;
          const int ifield = kv.iq;
          const auto& f3 = fields_3d(ie, ifield);
          const auto ef = [&] (const int& iedge, const int& k, const int& ip, const int& jp) {
            const auto& r3 = recv_3d_buffers(ie, ifield, iedge);
            // Using pointers here is a little bit of speedup that
            // can't be ignored in a kernel as important as un/pack.
            auto* const f3p = &f3(ip, jp, 0);
            const auto* const r3p = &r3(k, 0);
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
              [&] (const int& ilev) {
                f3p[ilev] += r3p[ilev];
              });
          };
          for (int k=0; k<NP; ++k) {
            ef(0, k, 0,    k);
            ef(1, k, NP-1, k);
            ef(2, k, k,    0);
            ef(3, k, k,    NP-1);
          }
          const auto cf = [&] (const int& icorner, const int& ip, const int& jp) {
            const auto& r3 = recv_3d_buffers(ie, ifield, icorner);
            if (r3.size() == 0)
              return;
            auto* const f3p = &f3(ip, jp, 0);
            const auto* const r3p = &r3(0, 0);
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
              [&] (const int& ilev) {
                f3p[ilev] += r3p[ilev];
              });
          };
          cf(4, 0,    0);
          cf(5, 0,    NP-1);
          cf(6, NP-1, 0);
          cf(7, NP-1, NP-1);
          if (rspheremp) {
            for (int i = 0; i < NP; ++i)
              for (int j = 0; j < NP; ++j) {
                auto* const f3p = &f3(i, j, 0);
                const auto& rsmp = (*rspheremp)(ie, i, j);
                Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
                  [&] (const int& ilev) {
                    f3p[ilev] *= rsmp;
                  });
              }
          }
        });
    }
  }
  // ...then unpack 3d interface fields (if any).
  if (m_num_3d_int_fields>0) {
    auto fields = m_3d_int_fields;
    auto recv_buffers = m_recv_3d_int_buffers;
    const auto num_fields = m_num_3d_int_fields;
    if (OnGpu<ExecSpace>::value) {
      const ConnectionHelpers helpers;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*num_fields*NUM_LEV_P),
        KOKKOS_LAMBDA(const int it) {
          const int ie = it / (num_fields*NUM_LEV_P);
          const int ifield = (it / NUM_LEV_P) % num_fields;
          const int ilev = it % NUM_LEV_P;
          const auto& f = fields(ie, ifield);
          for (int k=0; k<NP; ++k) {
            for (int iedge : helpers.UNPACK_EDGES_ORDER) {
              f(helpers.CONNECTION_PTS_FWD[iedge][k].ip,
                helpers.CONNECTION_PTS_FWD[iedge][k].jp, ilev)
                += recv_buffers(ie, ifield, iedge)(k, ilev);
            }
          }
          for (int icorner : helpers.UNPACK_CORNERS_ORDER) {
            if (recv_buffers(ie, ifield, icorner).size() > 0)
              f(helpers.CONNECTION_PTS_FWD[icorner][0].ip,
                helpers.CONNECTION_PTS_FWD[icorner][0].jp, ilev)
                += recv_buffers(ie, ifield, icorner)(0, ilev);
          }
        });
      if (rspheremp) {
        const auto rsmp = *rspheremp;
        Kokkos::parallel_for(
          Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*num_fields*NP*NP*NUM_LEV_P),
          KOKKOS_LAMBDA(const int it) {
            const int ie = it / (num_fields*NUM_LEV_P*NP*NP);
            const int ifield = (it / (NP*NP*NUM_LEV_P)) % num_fields;
            const int i = (it / (NP*NUM_LEV_P)) % NP;
            const int j = (it / NUM_LEV_P) % NP;
            const int ilev = it % NUM_LEV_P;
            fields(ie, ifield)(i, j, ilev) *= rsmp(ie, i, j);
          });
      }
    } else {
      const auto num_parallel_iterations = m_num_elems*num_fields;
      Kokkos::parallel_for(
        Kokkos::TeamPolicy<ExecSpace>(num_parallel_iterations, 1, NUM_LEV_P),
        KOKKOS_LAMBDA(const TeamMember& team) {
          Homme::KernelVariables kv(team, num_fields);
          const int ie = kv.ie;
          const int ifield = kv.iq;
          const auto& f = fields(ie, ifield);
          const auto ef = [&] (const int& iedge, const int& k, const int& ip, const int& jp) {
            const auto& rb = recv_buffers(ie, ifield, iedge);
            // Using pointers here is a little bit of speedup that
            // can't be ignored in a kernel as important as un/pack.
            auto* const fp = &f(ip, jp, 0);
            const auto* const rp = &rb(k, 0);
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV_P),
              [&] (const int& ilev) {
                fp[ilev] += rp[ilev];
              });
          };
          for (int k=0; k<NP; ++k) {
            ef(0, k, 0,    k);
            ef(1, k, NP-1, k);
            ef(2, k, k,    0);
            ef(3, k, k,    NP-1);
          }
          const auto cf = [&] (const int& icorner, const int& ip, const int& jp) {
            const auto& rb = recv_buffers(ie, ifield, icorner);
            if (rb.size() == 0) {
              return;
            }
            auto* const fp = &f(ip, jp, 0);
            const auto* const rp = &rb(0, 0);
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV_P),
              [&] (const int& ilev) {
                fp[ilev] += rp[ilev];
              });
          };
          cf(4, 0,    0);
          cf(5, 0,    NP-1);
          cf(6, NP-1, 0);
          cf(7, NP-1, NP-1);
          if (rspheremp) {
            for (int i = 0; i < NP; ++i)
              for (int j = 0; j < NP; ++j) {
                auto* const fp = &f(i, j, 0);
                const auto& rsmp = (*rspheremp)(ie, i, j);
                Kokkos::parallel_for(
                  Kokkos::ThreadVectorRange(kv.team, NUM_LEV_P),
                  [&] (const int& ilev) {
                    fp[ilev] *= rsmp;
                  });
              }
          }
        });
    }
  }
  ExecSpace::impl_static_fence();

  // If another BE structure starts an exchange, it has no way to check that
  // this object has finished its send requests, and may erroneously reuse the
  // buffers. Therefore, we must ensure that, upon return, all buffers are
  // reusable.

  tstart("be waitall 2");
  if ( ! m_send_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Waitall(m_send_requests.size(), m_send_requests.data(),
                                        MPI_STATUSES_IGNORE),
                            m_connectivity->get_comm().mpi_comm()); // Wait for all data to arrive
  tstop("be waitall 2");

  tstart("be recv_and_unpack book");
  // Release the send/recv buffers
  m_buffers_manager->unlock_buffers();
  m_send_pending = false;
  m_recv_pending = false;
  tstop("be recv_and_unpack book");
  tstop("be recv_and_unpack");
}

void BoundaryExchange::pack_and_send_min_max ()
{
  // The registration MUST be completed by now
  // Note: this also implies connectivity and buffers manager are valid
  assert (m_registration_completed);

  // Check that this object is setup to perform exchange_min_max and not exchange
  assert (m_exchange_type==MPI_EXCHANGE_MIN_MAX);

  // I am not sure why and if we could have this scenario, but just in case. I
  // think MPI *may* go bananas in this case
  if (m_num_1d_fields==0) {
    return;
  }

  // Check that buffers are not locked by someone else, then lock them
  assert (!m_buffers_manager->are_buffers_busy());
  m_buffers_manager->lock_buffers();

  // Check that this object is setup to perform exchange_min_max and not exchange
  assert (m_exchange_type==MPI_EXCHANGE_MIN_MAX);

  // If this is the first time we call this method, or if the MpiBuffersManager has performed a reallocation
  // since the last time this method was called, AND we are calling this method manually, without relying
  // on the exchange_min_max method to call it, then we need to rebuild all our internal buffer views
  if (!m_buffer_views_and_requests_built) {
    tstart("be build_buffer_views_and_requests");
    build_buffer_views_and_requests();
    tstop("be build_buffer_views_and_requests");
  }

  // NOTE: all of these temporary copies are necessary because of the issue of lambda function not
  //       capturing the this pointer correctly on the device.
  auto connections = m_connectivity->get_connections<ExecMemSpace>();
  auto fields_1d   = m_1d_fields;
  auto send_1d_buffers = m_send_1d_buffers;

  // ---- Pack ---- //
  const auto num_1d_fields = m_num_1d_fields;
  if (OnGpu<ExecSpace>::value) {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*m_num_1d_fields*NUM_CONNECTIONS*NUM_LEV),
      KOKKOS_LAMBDA(const int it) {
        const int ie = it / (num_1d_fields*NUM_CONNECTIONS*NUM_LEV);
        const int ifield = (it / (NUM_CONNECTIONS*NUM_LEV)) % num_1d_fields;
        const int iconn = (it / NUM_LEV) % NUM_CONNECTIONS;
        const int ilev = it % NUM_LEV;
        const ConnectionInfo& info = connections(ie, iconn);
        const LidGidPos& field_lidpos  = info.local;
        // For the buffer, in case of local connection, use remote info. In fact,
        // while with shared connections the mpi call will take care of "copying"
        // data to the remote recv buffer in the correct remote element lid, for
        // local connections we need to manually copy on the remote element
        // lid. We can do it here
        const LidGidPos& buffer_lidpos = info.sharing==etoi(ConnectionSharing::LOCAL) ? info.remote : info.local;

        send_1d_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos)(MAX_ID, ilev) =
          fields_1d(field_lidpos.lid, ifield)(MAX_ID, ilev);
        send_1d_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos)(MIN_ID, ilev) =
          fields_1d(field_lidpos.lid, ifield)(MIN_ID, ilev);
      });
  } else {
    const auto num_parallel_iterations = m_num_elems*num_1d_fields*NUM_CONNECTIONS;
    ThreadPreferences tp;
    tp.max_threads_usable = 1;
    tp.max_vectors_usable = NUM_LEV;
    const auto threads_vectors =
      DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
        num_parallel_iterations, tp);
    const auto policy = Kokkos::TeamPolicy<ExecSpace>(
      num_parallel_iterations, threads_vectors.first, threads_vectors.second);
    Kokkos::parallel_for(
      policy,
      KOKKOS_LAMBDA(const TeamMember& team) {
        Homme::KernelVariables kv(team, num_1d_fields*NUM_CONNECTIONS);
        const int ie = kv.ie;
        const int iconn = kv.iq / num_1d_fields;
        const int ifield = kv.iq % num_1d_fields;
        const ConnectionInfo& info = connections(ie, iconn);
        if (info.kind == etoi(ConnectionSharing::MISSING)) return;
        const LidGidPos& field_lidpos = info.local;
        const LidGidPos& buffer_lidpos = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                          info.remote :
                                          info.local);
        const auto& sb = send_1d_buffers(buffer_lidpos.lid, ifield, buffer_lidpos.pos);
        const auto& f1 = fields_1d(field_lidpos.lid, ifield);
        for (int k = 0; k < 2; ++k) {
          auto* const sbp = &sb(k, 0);
          const auto* const f1p = &f1(k, 0);
          Kokkos::parallel_for(
            Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
            [&] (const int& ilev) {
              sbp[ilev] = f1p[ilev];
            });
        }
      });
  }
  ExecSpace::impl_static_fence();

  // ---- Send ---- //
  m_buffers_manager->sync_send_buffer(this);
  if ( ! m_send_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Startall(m_send_requests.size(), m_send_requests.data()),
                            m_connectivity->get_comm().mpi_comm());

  // Mark send buffer as busy
  m_send_pending = true;
}

void BoundaryExchange::recv_and_unpack_min_max ()
{
  // The registration MUST be completed by now
  // Note: this also implies connectivity and buffers manager are valid
  assert (m_registration_completed);

  // Check that this object is setup to perform exchange_min_max and not exchange
  assert (m_exchange_type==MPI_EXCHANGE_MIN_MAX);

  // I am not sure why and if we could have this scenario, but just in case. I
  // think MPI *may* go bananas in this case
  if (m_num_1d_fields==0) {
    return;
  }

  // If I am doing pack_and_send and recv_and_unpack manually (rather than through 'exchange'),
  // then I need to start receiving now (otherwise it is done already inside 'exchange')
  if (!m_recv_pending) {
    // If you are doing send/recv manually, don't call recv without a send, or
    // else you'll be stuck waiting later on
    assert (m_send_pending);

    if ( ! m_recv_requests.empty())
      HOMMEXX_MPI_CHECK_ERROR(MPI_Startall(m_recv_requests.size(), m_recv_requests.data()),
                              m_connectivity->get_comm().mpi_comm());
    m_recv_pending = true;
  }

  // ---- Recv ---- //
  if ( ! m_recv_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Waitall(m_recv_requests.size(), m_recv_requests.data(), MPI_STATUSES_IGNORE),
                            m_connectivity->get_comm().mpi_comm()); // Wait for all data to arrive

  m_buffers_manager->sync_recv_buffer(this); // Deep copy mpi_recv_buffer into recv_buffer (no op if MPI is on device)

  // NOTE: all of these temporary copies are necessary because of the issue of lambda function not
  //       capturing the this pointer correctly on the device.
  auto connections = m_connectivity->get_connections<ExecMemSpace>();
  auto fields_1d = m_1d_fields;
  auto recv_1d_buffers = m_recv_1d_buffers;

  // --- Unpack --- //
  const auto num_1d_fields = m_num_1d_fields;
  if (OnGpu<ExecSpace>::value) {
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*m_num_1d_fields*NUM_LEV),
      KOKKOS_LAMBDA(const int it) {
        const int ie = it / (num_1d_fields*NUM_LEV);
        const int ifield = (it / NUM_LEV) % num_1d_fields;
        const int ilev = it % NUM_LEV;
        for (int neighbor=0; neighbor<NUM_CONNECTIONS; ++neighbor) {
          // Note: for min/max exchange, we really need to skip MISSING
          //       connections (while for 'normal' exchange, the missing recv
          //       buffer points to a blackhole fileld with 0's, which do not
          //       alter the accummulation)
          if (connections(ie, neighbor).kind==etoi(ConnectionKind::MISSING)) {
            continue;
          }
          fields_1d(ie, ifield)(MAX_ID, ilev) = max(fields_1d(ie, ifield)(MAX_ID, ilev),
                                                    recv_1d_buffers(ie, ifield, neighbor)(MAX_ID, ilev));
          fields_1d(ie, ifield)(MIN_ID, ilev) = min(fields_1d(ie, ifield)(MIN_ID, ilev),
                                                    recv_1d_buffers(ie, ifield, neighbor)(MIN_ID, ilev));
        }
      });
  } else {
    const auto num_parallel_iterations = m_num_elems*num_1d_fields;
    ThreadPreferences tp;
    tp.max_threads_usable = 1;
    tp.max_vectors_usable = NUM_LEV;
    const auto threads_vectors =
      DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
        num_parallel_iterations, tp);
    const auto policy = Kokkos::TeamPolicy<ExecSpace>(
      num_parallel_iterations, threads_vectors.first, threads_vectors.second);
    Kokkos::parallel_for(
      policy,
      KOKKOS_LAMBDA(const TeamMember& team) {
        Homme::KernelVariables kv(team, num_1d_fields);
        const int ie = kv.ie;
        const int ifield = kv.iq;
        for (int iconn = 0; iconn < NUM_CONNECTIONS; ++iconn) {
          const ConnectionInfo& info = connections(ie, iconn);
          if (info.kind == etoi(ConnectionSharing::MISSING)) continue;
          const auto& rb = recv_1d_buffers(ie, ifield, iconn);
          const auto& f1 = fields_1d(ie, ifield);
          {
            const auto* const rbp = &rb(0, 0);
            auto* const f1p = &f1(0, 0);
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
              [&] (const int& ilev) {
                f1p[ilev] = min(f1p[ilev], rbp[ilev]);
              });
          }
          {
            const auto* const rbp = &rb(1, 0);
            auto* const f1p = &f1(1, 0);
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(kv.team, NUM_LEV),
              [&] (const int& ilev) {
                f1p[ilev] = max(f1p[ilev], rbp[ilev]);
              });
          }
        }
      });
  }
  ExecSpace::impl_static_fence();

  // If another BE structure starts an exchange, it has no way to check that
  // this object has finished its send requests, and may erroneously reuse the
  // buffers. Therefore, we must ensure that, upon return, all buffers are
  // reusable.
  if ( ! m_send_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Waitall(m_send_requests.size(), m_send_requests.data(), MPI_STATUSES_IGNORE),
                            m_connectivity->get_comm().mpi_comm()); // Wait for all data to arrive

  // Release the send/recv buffers
  m_buffers_manager->unlock_buffers();
  m_send_pending = false;
  m_recv_pending = false;
}

void BoundaryExchange
::free_requests () {
  for (size_t i=0; i<m_send_requests.size(); ++i)
    HOMMEXX_MPI_CHECK_ERROR(MPI_Request_free(&m_send_requests[i]),
                            m_connectivity->get_comm().mpi_comm());
  m_send_requests.clear();
  for (size_t i=0; i<m_recv_requests.size(); ++i)
    HOMMEXX_MPI_CHECK_ERROR(MPI_Request_free(&m_recv_requests[i]),
                            m_connectivity->get_comm().mpi_comm());
  m_recv_requests.clear();
}

// A slot is the space in a communication buffer for an (element, connection)
// pair. The slot index space numbers slots so that, first, they are contiguous
// by remote PID and, second, within a PID block, each comm partner agrees on
// order of slots.
void BoundaryExchange
::init_slot_idx_to_elem_conn_pair (
  std::vector<int>& slot_idx_to_elem_conn_pair,
  std::vector<int>& pids, std::vector<int>& pid_offsets)
{
  struct IP {
    int i, ord, pid;
    bool operator< (const IP& o) const {
      if (pid < o.pid) return true;
      if (pid > o.pid) return false;
      return ord < o.ord;
    }
  };
  std::vector<IP> i2remote(m_num_elems*NUM_CONNECTIONS);

  const auto& connections = m_connectivity->get_connections<HostMemSpace>();
  for (int ie = 0; ie < m_num_elems; ++ie)
    for (int iconn = 0; iconn < NUM_CONNECTIONS; ++iconn) {
      const auto& info = connections(ie, iconn);
      const int k = ie*NUM_CONNECTIONS + iconn;
      auto& i2r = i2remote[k];
      // Original sequence through (element, connection) pairs.
      i2r.i = k;
      // An ordering of the message buffer upon which both members of the
      // communication pair agree.
      if (info.local.gid < info.remote.gid)
        i2r.ord = info.local.gid*NUM_CONNECTIONS + info.local.pos;
      else
        i2r.ord = info.remote.gid*NUM_CONNECTIONS + info.remote.pos;
      // If local, indicate with -1, which is < the smallest pid of 0.
      i2r.pid = -1;
      if (info.sharing != etoi(ConnectionSharing::SHARED)) continue;
      i2r.pid = info.remote_pid;
    }

  // Sort so that, first, all (element, connection) pairs having the same
  // remote_pid are contiguous; second, within such a block, ord is
  // ascending. The first lets us set up comm buffers so monolithic messages
  // slot right in. The second means that the send and recv partners agree on
  // how the monolithic message is packed.
  std::sort(i2remote.begin(), i2remote.end());

  // Collect the unique remote_pids and get the offsets of the contiguous blocks
  // of them.
  slot_idx_to_elem_conn_pair.resize(m_num_elems*NUM_CONNECTIONS);
  pids.clear();
  pid_offsets.clear();
  int prev_pid = -2;
  for (int k = 0; k < m_num_elems*NUM_CONNECTIONS; ++k) {
    const auto& i2r = i2remote[k];
    if (i2r.pid > prev_pid && i2r.pid != -1) {
      pids.push_back(i2r.pid);
      pid_offsets.push_back(k);
      prev_pid = i2r.pid;
    }
    const int ie = i2r.i / NUM_CONNECTIONS;
    const int iconn = i2r.i % NUM_CONNECTIONS;
    slot_idx_to_elem_conn_pair[k] = ie*NUM_CONNECTIONS + iconn;
  }
  pid_offsets.push_back(m_num_elems*NUM_CONNECTIONS);
}

void BoundaryExchange::clear_buffer_views_and_requests ()
{
  // MpiBuffersManager calls this method upon (re)allocation of buffers, so that all its customers are forced to
  // recompute their internal buffers views. However, if the views were not yet built, we can skip this
  if (!m_buffer_views_and_requests_built) {
    return;
  }

  // The connectivity must be valid here
  assert (m_connectivity);

  // Destroy each request
  free_requests();

  // Clear buffer views
  m_send_1d_buffers = decltype(m_send_1d_buffers)("m_send_1d_buffers", 0, 0);
  m_recv_1d_buffers = decltype(m_recv_1d_buffers)("m_recv_1d_buffers", 0, 0);
  m_send_2d_buffers = decltype(m_send_2d_buffers)("m_send_2d_buffers", 0, 0);
  m_recv_2d_buffers = decltype(m_recv_2d_buffers)("m_recv_2d_buffers", 0, 0);
  m_send_3d_buffers = decltype(m_send_3d_buffers)("m_send_3d_buffers", 0, 0);
  m_recv_3d_buffers = decltype(m_recv_3d_buffers)("m_recv_3d_buffers", 0, 0);
  m_send_3d_int_buffers = decltype(m_send_3d_int_buffers)("m_send_3d_int_buffers", 0, 0);
  m_recv_3d_int_buffers = decltype(m_recv_3d_int_buffers)("m_recv_3d_int_buffers", 0, 0);

  // Done
  m_buffer_views_and_requests_built = false;
}

void BoundaryExchange::waitall()
{
  if (!m_send_pending && !m_recv_pending) {
    return;
  }

  // At this point, the connectivity MUST be valid
  assert (m_connectivity);

  // Safety check
  assert (m_buffers_manager->are_buffers_busy());

  if ( ! m_send_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Waitall(m_send_requests.size(), m_send_requests.data(), MPI_STATUSES_IGNORE),
                            m_connectivity->get_comm().mpi_comm());
  if ( ! m_recv_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Waitall(m_recv_requests.size(), m_recv_requests.data(), MPI_STATUSES_IGNORE),
                            m_connectivity->get_comm().mpi_comm());

  m_buffers_manager->unlock_buffers();
}

} // namespace Homme
