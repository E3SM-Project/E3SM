/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "BoundaryExchange.hpp"

#include "BuffersManager.hpp"
#include "KernelVariables.hpp"
#include "profiling.hpp"

#include "utilities/VectorUtils.hpp"

#define tstart(x)
#define tstop(x)

namespace Homme
{

// ======================== IMPLEMENTATION ======================== //

BoundaryExchange::BoundaryExchange()
{
  m_num_1d_fields = 0;
  m_num_2d_fields = 0;
  m_num_3d_fields = 0;

  m_connectivity    = std::shared_ptr<Connectivity>();
  m_buffers_manager = std::shared_ptr<BuffersManager>();

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

BoundaryExchange::BoundaryExchange(std::shared_ptr<Connectivity> connectivity, std::shared_ptr<BuffersManager> buffers_manager)
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

void BoundaryExchange::set_buffers_manager (std::shared_ptr<BuffersManager> buffers_manager)
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

void BoundaryExchange::set_num_fields (const int num_1d_fields, const int num_2d_fields, const int num_3d_fields)
{
  // We don't allow to call this method twice in a row. If you want to change the number of fields,
  // you need to call clean_up first, to get a fresh new BoundaryExchange.
  // Note: if you do call clean_up and then again set_num_field and the new number of fields
  //       are smaller than the previous ones, BuffersManager will NOT shrink the buffers, so
  //       you may be left with buffers that are larger than what you need.
  assert (m_cleaned_up);

  // Make sure the connectivity is valid: must at least be a valid pointer, and be initialized, i.e.,
  // store a valid number of elements, but may be finalized later (before registration_completed call though)
  assert (m_connectivity && m_connectivity->is_initialized());

  // We strongly advocate for not using the same BE object for both 'standard' exchange and min/max exchange
  assert (!(num_1d_fields>0 && (num_2d_fields>0 || num_2d_fields>0)));

  // Note: we do not set m_num_1d_fields, m_num_2d_fields and m_num_3d_fields, since we will use them as
  //       progressive indices while adding fields. Then, during registration_completed,
  //       we will check that they match the 2nd dimension of m_1d_fields, m_2d_fields and m_3d_fields.

  // Create the fields views
  m_1d_fields = decltype(m_1d_fields)("1d fields", m_num_elems, num_1d_fields);
  m_2d_fields = decltype(m_2d_fields)("2d fields", m_num_elems, num_2d_fields);
  m_3d_fields = decltype(m_3d_fields)("3d fields", m_num_elems, num_3d_fields);

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
  m_2d_fields = decltype(m_2d_fields)("m_2d_fields", 0, 0);
  m_3d_fields = decltype(m_3d_fields)("m_3d_fields", 0, 0);

  m_num_2d_fields = 0;
  m_num_3d_fields = 0;

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
  m_elem_buf_size[etoi(ConnectionKind::CORNER)] = m_num_1d_fields*2*NUM_LEV*VECTOR_SIZE + (m_num_2d_fields + m_num_3d_fields*NUM_LEV*VECTOR_SIZE) * 1;
  m_elem_buf_size[etoi(ConnectionKind::EDGE)]   = m_num_1d_fields*2*NUM_LEV*VECTOR_SIZE + (m_num_2d_fields + m_num_3d_fields*NUM_LEV*VECTOR_SIZE) * NP;

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
  if (m_num_2d_fields+m_num_3d_fields==0) {
    return;
  }

  // If this is the first time we call the exchange method, or if the BuffersManager has performed a reallocation
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

  // If this is the first time we call the exchange method, or if the BuffersManager has performed a reallocation
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
  if (m_num_2d_fields+m_num_3d_fields==0) {
    return;
  }

  // Check that buffers are not locked by someone else, then lock them
  assert (!m_buffers_manager->are_buffers_busy());
  m_buffers_manager->lock_buffers();

  // If this is the first time we call this method, or if the BuffersManager has performed a reallocation
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
  // ...then pack 3d fields (if any).
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
  ExecSpace::fence();

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
  // ...then unpack 3d fields.
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
        const auto r = *rspheremp;
        Kokkos::parallel_for(
          Kokkos::RangePolicy<ExecSpace>(0, m_num_elems*m_num_3d_fields*NP*NP*NUM_LEV),
          KOKKOS_LAMBDA(const int it) {
            const int ie = it / (num_3d_fields*NUM_LEV*NP*NP);
            const int ifield = (it / (NP*NP*NUM_LEV)) % num_3d_fields;
            const int i = (it / (NP*NUM_LEV)) % NP;
            const int j = (it / NUM_LEV) % NP;
            const int ilev = it % NUM_LEV;
            fields_3d(ie, ifield)(i, j, ilev) *= r(ie, i, j);
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
  ExecSpace::fence();

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

  // If this is the first time we call this method, or if the BuffersManager has performed a reallocation
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
  ExecSpace::fence();

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
  ExecSpace::fence();

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

void BoundaryExchange::build_buffer_views_and_requests()
{
  // If we already set the buffers before, then nothing to be done here
  if (m_buffer_views_and_requests_built) {
    return;
  }

  // Check that the BuffersManager is present and was setup with enough storage
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
    }
  }
  Kokkos::deep_copy(m_send_1d_buffers, h_send_1d_buffers);
  Kokkos::deep_copy(m_send_2d_buffers, h_send_2d_buffers);
  Kokkos::deep_copy(m_send_3d_buffers, h_send_3d_buffers);
  Kokkos::deep_copy(m_recv_1d_buffers, h_recv_1d_buffers);
  Kokkos::deep_copy(m_recv_2d_buffers, h_recv_2d_buffers);
  Kokkos::deep_copy(m_recv_3d_buffers, h_recv_3d_buffers);

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
  // BuffersManager calls this method upon (re)allocation of buffers, so that all its customers are forced to
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
