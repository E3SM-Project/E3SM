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

  int single_ptr_buf_size = m_num_2d_fields + m_num_3d_int_fields*NUM_LEV_P*VECTOR_SIZE;
  for (int i = 0; i < m_num_3d_fields; ++i)
    single_ptr_buf_size += m_3d_nlev_pack[i]*VECTOR_SIZE;
  m_elem_buf_size[etoi(ConnectionKind::CORNER)] = m_num_1d_fields*2*NUM_LEV*VECTOR_SIZE + single_ptr_buf_size * 1;
  m_elem_buf_size[etoi(ConnectionKind::EDGE)]   = m_num_1d_fields*2*NUM_LEV*VECTOR_SIZE + single_ptr_buf_size * NP;

  // Determine what kind of BE is this (exchange or exchange_min_max)
  m_exchange_type = m_num_1d_fields>0 ? MPI_EXCHANGE_MIN_MAX : MPI_EXCHANGE;

  // Finalize bookkeeping for any exchange on fewer than NUM_LEV levels.
  {
    bool need_nlev_pack = false;
    for (int i = 0; i < m_num_3d_fields; ++i)
      if (m_3d_nlev_pack[i] != NUM_LEV) {
        Errors::runtime_check(m_3d_nlev_pack[i] < NUM_LEV,
                              "Optional nlev must be <= NUM_LEV");
        Errors::runtime_check(m_3d_nlev_pack[i] > 0,
                              "Optional nlev must be > 0");
        need_nlev_pack = true;
        break;
      }
    if (need_nlev_pack) {
      m_3d_nlev_pack_d = ExecViewManaged<int*>("m_3d_nlev_pack_d", m_3d_nlev_pack.size());
      const auto h = Kokkos::create_mirror_view(m_3d_nlev_pack_d);
      for (int i = 0; i < m_num_3d_fields; ++i) h(i) = m_3d_nlev_pack[i];
      Kokkos::deep_copy(m_3d_nlev_pack_d, h);
    }
    // Clear the host vector if it's not needed.
    if ( ! need_nlev_pack) m_3d_nlev_pack = decltype(m_3d_nlev_pack)();
  }

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

static void
pack (const ExecViewUnmanaged<const HaloExchangeUnstructuredConnectionInfo*> ucon,
      const ExecViewUnmanaged<const int*> ucon_ptr,
      const ExecViewUnmanaged<ExecViewManaged<Real[NP][NP]>**> fields_2d,
      const ExecViewUnmanaged<ExecViewUnmanaged<Real*>**> send_2d_buffers,
      const int num_elems, const int num_2d_fields) {
  HOMMEXX_STATIC const ConnectionHelpers helpers;
  const int nconn = ucon.extent_int(0);
  Kokkos::parallel_for(
    Kokkos::RangePolicy<ExecSpace>(0, num_2d_fields*nconn),
    KOKKOS_LAMBDA(const int it) {
      const int iconn = it / num_2d_fields;
      const int ifield = it % num_2d_fields;
      const auto& info = ucon(iconn);
      const int buffer_iconn = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                info.sharing_local_remote_iconn :
                                iconn);
      const auto& pts = helpers.CONNECTION_PTS[info.direction][info.local_dir];
      const auto& sb = send_2d_buffers(ifield, buffer_iconn);
      const auto& f2 = fields_2d(info.local_lid, ifield);
      for (int k = 0; k < helpers.CONNECTION_SIZE[info.kind]; ++k)
        sb(k) = f2(pts[k].ip, pts[k].jp);
    });
}

template <int NUM_LEV_PACKS, bool partial_column=false>
static void
pack (const ExecViewUnmanaged<const HaloExchangeUnstructuredConnectionInfo*> ucon,
      const ExecViewUnmanaged<const int*> ucon_ptr,
      const ExecViewUnmanaged<ExecViewManaged<Scalar[NP][NP][NUM_LEV_PACKS]>**> fields_3d,
      const ExecViewUnmanaged<ExecViewUnmanaged<Scalar**>**> send_3d_buffers,
      const int num_elems, const int num_3d_fields,
      ExecViewManaged<int*>* nlev_packs_ = nullptr) {
  assert(partial_column == (nlev_packs_ != nullptr));
  if (partial_column) assert(nlev_packs_->extent_int(0) == num_3d_fields);
  ExecViewUnmanaged<const int*> nlev_packs;
  if (partial_column) nlev_packs = *nlev_packs_;
  if (OnGpu<ExecSpace>::value) {
    const ConnectionHelpers helpers;
    const int nconn = ucon.extent_int(0);
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, num_3d_fields*nconn*NUM_LEV_PACKS),
      KOKKOS_LAMBDA(const int it) {
        const int ilev = it % NUM_LEV_PACKS;
        const int ifield = (it / NUM_LEV_PACKS) % num_3d_fields;
        if (partial_column) { // compile out if !partial_column
          if (ilev >= nlev_packs(ifield))
            return;
        }
        const int iconn = it / (num_3d_fields*NUM_LEV_PACKS);
        const auto& info = ucon(iconn);
        const int buffer_iconn = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                  info.sharing_local_remote_iconn :
                                  iconn);
        const auto& pts = helpers.CONNECTION_PTS[info.direction][info.local_dir];
        const auto& sb = send_3d_buffers(ifield, buffer_iconn);
        const auto& f3 = fields_3d(info.local_lid, ifield);
        for (int k = 0; k < helpers.CONNECTION_SIZE[info.kind]; ++k)
          sb(k, ilev) = f3(pts[k].ip, pts[k].jp, ilev);
      });
  } else {
    const auto num_parallel_iterations = num_elems*num_3d_fields;
    ThreadPreferences tp;
    tp.max_threads_usable = NP;
    tp.max_vectors_usable = NUM_LEV_PACKS;
    const auto threads_vectors =
      DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
        num_parallel_iterations, tp);
    const auto policy = Kokkos::TeamPolicy<ExecSpace>(
      num_parallel_iterations, threads_vectors.first, threads_vectors.second);
    HOMMEXX_STATIC const ConnectionHelpers helpers;
    Kokkos::parallel_for(policy,
      KOKKOS_LAMBDA(const TeamMember& team) {
        Homme::KernelVariables kv(team, num_3d_fields);
        const int ie = kv.ie;
        const int ifield = kv.iq;
        const auto tvr = Kokkos::ThreadVectorRange(
          kv.team, partial_column ? nlev_packs(ifield) : NUM_LEV_PACKS);
        const int iconn_end = ucon_ptr(ie+1);
        for (int iconn = ucon_ptr(ie); iconn < iconn_end; ++iconn) {
          const auto& info = ucon(iconn);
          assert(info.kind != etoi(ConnectionSharing::MISSING));
          const int buffer_iconn = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                    info.sharing_local_remote_iconn :
                                    iconn);
          const auto& pts = helpers.CONNECTION_PTS[info.direction][info.local_dir];
          const auto& sb = send_3d_buffers(ifield, buffer_iconn);
          assert(info.local_lid == ie);
          const auto& f3 = fields_3d(ie, ifield);
          Kokkos::parallel_for(
            Kokkos::TeamThreadRange(kv.team, helpers.CONNECTION_SIZE[info.kind]),
            [&] (const int& k) {
              auto* const sbp = &sb(k, 0);
              const auto* const f3p = &f3(pts[k].ip, pts[k].jp, 0);
              Kokkos::parallel_for(tvr, [&] (const int& ilev) { sbp[ilev] = f3p[ilev]; });
            });
        }
      });
  }
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
  const auto& ucon = m_connectivity->get_d_ucon();
  const auto& ucon_ptr = m_connectivity->get_d_ucon_ptr();
  // First, pack 2d fields (if any)...
  if (m_num_2d_fields > 0)
    pack(ucon, ucon_ptr, m_2d_fields, m_send_2d_buffers, m_num_elems,
         m_num_2d_fields);
  // ...then pack 3d fields (if any)...
  if (m_num_3d_fields > 0) {
    if (m_3d_nlev_pack_d.size() > 0)
      pack<NUM_LEV, true>(ucon, ucon_ptr, m_3d_fields, m_send_3d_buffers,
                          m_num_elems, m_num_3d_fields, &m_3d_nlev_pack_d);
    else
      pack<NUM_LEV>(ucon, ucon_ptr, m_3d_fields, m_send_3d_buffers,
                    m_num_elems, m_num_3d_fields);
  }
  // ...then pack 3d interface fields (if any)
  if (m_num_3d_int_fields > 0)
    pack<NUM_LEV_P>(ucon, ucon_ptr, m_3d_int_fields, m_send_3d_int_buffers,
                    m_num_elems, m_num_3d_int_fields);
  Kokkos::fence();

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

// assume:conn-edges-snwe
static void
unpack (const ExecViewUnmanaged<const HaloExchangeUnstructuredConnectionInfo*> ucon,
        const ExecViewUnmanaged<const int*> ucon_ptr,
        const ExecViewUnmanaged<ExecViewManaged<Real[NP][NP]>**> fields_2d,
        const ExecViewUnmanaged<ExecViewUnmanaged<Real*>**> recv_2d_buffers,
        const ExecViewUnmanaged<const Real * [NP][NP]>* rspheremp,
        const int num_elems, const int num_2d_fields) {
  HOMMEXX_STATIC const ConnectionHelpers helpers;
  Kokkos::parallel_for(
    Kokkos::RangePolicy<ExecSpace>(0, num_elems*num_2d_fields),
    KOKKOS_LAMBDA(const int it) {
      const int ie = it / num_2d_fields;
      const int ifield = it % num_2d_fields;
      const auto iconn_beg = ucon_ptr(ie), iconn_end = ucon_ptr(ie+1);
      const auto& f2 = fields_2d(ie, ifield);
      for (int k = 0; k < NP; ++k) {
        for (const int iedge : helpers.UNPACK_EDGES_ORDER) {
          f2(helpers.CONNECTION_PTS_FWD[iedge][k].ip,
             helpers.CONNECTION_PTS_FWD[iedge][k].jp)
            += recv_2d_buffers(ifield, iconn_beg + iedge)(k);
        }
      }
      for (int iconn = iconn_beg + 4; iconn < iconn_end; ++iconn) {
        const auto dir = ucon(iconn).local_dir;
        f2(helpers.CONNECTION_PTS_FWD[dir][0].ip,
           helpers.CONNECTION_PTS_FWD[dir][0].jp)
          += recv_2d_buffers(ifield, iconn)(0);
      }
    });  
  if (rspheremp) {
    Kokkos::fence();
    const auto rsmp = *rspheremp;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, num_elems*num_2d_fields*NP*NP),
      KOKKOS_LAMBDA(const int it) {
        const int ie = it / (num_2d_fields*NP*NP);
        const int ifield = (it / (NP*NP)) % num_2d_fields;
        const int i = (it / NP) % NP;
        const int j = it % NP;
        fields_2d(ie, ifield)(i, j) *= rsmp(ie, i, j);
      });
  }
}

// assume:conn-edges-snwe
template <int NUM_LEV_PACKS, bool partial_column=false>
static void
unpack (const ExecViewUnmanaged<const HaloExchangeUnstructuredConnectionInfo*> ucon,
        const ExecViewUnmanaged<const int*> ucon_ptr,
        const ExecViewUnmanaged<ExecViewManaged<Scalar[NP][NP][NUM_LEV_PACKS]>**> fields_3d,
        const ExecViewUnmanaged<ExecViewUnmanaged<Scalar**>**> recv_3d_buffers,
        const ExecViewUnmanaged<const Real * [NP][NP]>* rspheremp,
        const int num_elems, const int num_3d_fields,
        ExecViewManaged<int*>* nlev_packs_ = nullptr) {
  assert(partial_column == (nlev_packs_ != nullptr));
  if (partial_column) assert(nlev_packs_->extent_int(0) == num_3d_fields);
  ExecViewUnmanaged<const int*> nlev_packs;
  if (partial_column) nlev_packs = *nlev_packs_;
  if (OnGpu<ExecSpace>::value) {
    const ConnectionHelpers helpers;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, num_elems*num_3d_fields*NUM_LEV_PACKS),
      KOKKOS_LAMBDA(const int it) {
        const int ifield = (it / NUM_LEV_PACKS) % num_3d_fields;
        const int ilev = it % NUM_LEV_PACKS;
        if (partial_column) { // compile out if !partial_column
          if (ilev >= nlev_packs(ifield))
            return;
        }
        const int ie = it / (num_3d_fields*NUM_LEV_PACKS);
        const auto iconn_beg = ucon_ptr(ie);
        const auto& f3 = fields_3d(ie, ifield);
        for (int k = 0; k < NP; ++k) {
          for (const int iedge : helpers.UNPACK_EDGES_ORDER) {
            const auto& pts = helpers.CONNECTION_PTS_FWD[iedge][k];
            f3(pts.ip, pts.jp, ilev) +=
              recv_3d_buffers(ifield, iconn_beg + iedge)(k, ilev);
          }
        }
        const auto iconn_end = ucon_ptr(ie+1);
        for (int iconn = iconn_beg + 4; iconn < iconn_end; ++iconn) {
          const auto& pts = helpers.CONNECTION_PTS_FWD[ucon(iconn).local_dir][0];
          f3(pts.ip, pts.jp, ilev) +=
            recv_3d_buffers(ifield, iconn)(0, ilev);
        }
      });
    if (rspheremp) {
      Kokkos::fence();
      const auto rsmp = *rspheremp;
      Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, num_elems*num_3d_fields*NP*NP*NUM_LEV_PACKS),
        KOKKOS_LAMBDA(const int it) {
          const int ie = it / (num_3d_fields*NUM_LEV_PACKS*NP*NP);
          const int ifield = (it / (NP*NP*NUM_LEV_PACKS)) % num_3d_fields;
          const int i = (it / (NP*NUM_LEV_PACKS)) % NP;
          const int j = (it / NUM_LEV_PACKS) % NP;
          const int ilev = it % NUM_LEV_PACKS;
          fields_3d(ie, ifield)(i, j, ilev) *= rsmp(ie, i, j);
        });
    }
  } else {
    HOMMEXX_STATIC const ConnectionHelpers helpers;
    const auto num_parallel_iterations = num_elems*num_3d_fields;
    Kokkos::parallel_for(
      Kokkos::TeamPolicy<ExecSpace>(num_parallel_iterations, 1, NUM_LEV_PACKS),
      KOKKOS_LAMBDA(const TeamMember& team) {
        Homme::KernelVariables kv(team, num_3d_fields);
        const int ie = kv.ie;
        const int ifield = kv.iq;
        const auto tvr = Kokkos::ThreadVectorRange(
          kv.team, partial_column ? nlev_packs(ifield) : NUM_LEV_PACKS);
        const auto& f3 = fields_3d(ie, ifield);
        const auto iconn_beg = ucon_ptr(ie), iconn_end = ucon_ptr(ie+1);
        const auto ef = [&] (const int& iedge, const int& k, const int& ip, const int& jp) {
          const auto& r3 = recv_3d_buffers(ifield, iconn_beg + iedge);
          auto* const f3p = &f3(ip, jp, 0);
          const auto* const r3p = &r3(k, 0);
          Kokkos::parallel_for(tvr, [&] (const int& ilev) { f3p[ilev] += r3p[ilev]; });
        };
        for (int k = 0; k < NP; ++k) {
          ef(0, k, 0,    k   );
          ef(1, k, NP-1, k   );
          ef(2, k, k,    0   );
          ef(3, k, k,    NP-1);
        }
        for (int iconn = iconn_beg + 4; iconn < iconn_end; ++iconn) {
          const auto dir = ucon(iconn).local_dir;
          const auto& r3 = recv_3d_buffers(ifield, iconn);
          auto* const f3p = &f3(helpers.CONNECTION_PTS_FWD[dir][0].ip,
                                helpers.CONNECTION_PTS_FWD[dir][0].jp, 0);
          assert(r3.size() > 0);
          const auto* const r3p = &r3(0, 0);
          Kokkos::parallel_for(tvr, [&] (const int& ilev) { f3p[ilev] += r3p[ilev]; });
        }
        if (rspheremp) {
          for (int i = 0; i < NP; ++i)
            for (int j = 0; j < NP; ++j) {
              auto* const f3p = &f3(i, j, 0);
              const auto& rsmp = (*rspheremp)(ie, i, j);
              Kokkos::parallel_for(tvr, [&] (const int& ilev) { f3p[ilev] *= rsmp; });
            }
        }
      });
  }
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
  const auto& ucon = m_connectivity->get_d_ucon();
  const auto& ucon_ptr = m_connectivity->get_d_ucon_ptr();
  // First, unpack 2d fields (if any)...
  if (m_num_2d_fields>0)
    unpack(ucon, ucon_ptr, m_2d_fields, m_recv_2d_buffers, rspheremp, m_num_elems,
           m_num_2d_fields);
  // ...then unpack 3d fields (if any)...
  if (m_num_3d_fields>0) {
    if (m_3d_nlev_pack_d.size() > 0)
      unpack<NUM_LEV, true>(ucon, ucon_ptr, m_3d_fields, m_recv_3d_buffers, rspheremp,
                            m_num_elems, m_num_3d_fields, &m_3d_nlev_pack_d);
    else
      unpack<NUM_LEV>(ucon, ucon_ptr, m_3d_fields, m_recv_3d_buffers, rspheremp,
                      m_num_elems, m_num_3d_fields);
  }
  // ...then unpack 3d interface fields (if any).
  if (m_num_3d_int_fields > 0)
    unpack<NUM_LEV_P>(ucon, ucon_ptr, m_3d_int_fields, m_recv_3d_int_buffers, rspheremp,
                      m_num_elems, m_num_3d_int_fields);
  Kokkos::fence();

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

static void pack_min_max (
  const ExecViewUnmanaged<const HaloExchangeUnstructuredConnectionInfo*> ucon,
  const ExecViewUnmanaged<const int*> ucon_ptr,
  const ExecViewUnmanaged<ExecViewManaged<Scalar[2][NUM_LEV]>**> fields_1d,
  const ExecViewUnmanaged<ExecViewUnmanaged<Scalar[2][NUM_LEV]>**> send_1d_buffers,
  const int num_elems, const int num_1d_fields)
{
  if (OnGpu<ExecSpace>::value) {
    const ConnectionHelpers helpers;
    const int nconn = ucon.extent_int(0);
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, num_1d_fields*nconn*NUM_LEV),
      KOKKOS_LAMBDA(const int it) {
        const int iconn = it / (num_1d_fields*NUM_LEV);
        const int ifield = (it / NUM_LEV) % num_1d_fields;
        const int ilev = it % NUM_LEV;
        const auto& info = ucon(iconn);
        const int buffer_iconn = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                  info.sharing_local_remote_iconn :
                                  iconn);
        const auto& sb = send_1d_buffers(ifield, buffer_iconn);
        const auto& f1 = fields_1d(info.local_lid, ifield);
        for (int k = 0; k < 2; ++k)
          sb(k, ilev) = f1(k, ilev);
      });
  } else {
    const auto num_parallel_iterations = num_elems*num_1d_fields;
    ThreadPreferences tp;
    tp.max_threads_usable = 2;
    tp.max_vectors_usable = NUM_LEV;
    const auto threads_vectors =
      DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
        num_parallel_iterations, tp);
    const auto policy = Kokkos::TeamPolicy<ExecSpace>(
      num_parallel_iterations, threads_vectors.first, threads_vectors.second);
    HOMMEXX_STATIC const ConnectionHelpers helpers;
    Kokkos::parallel_for(policy,
      KOKKOS_LAMBDA(const TeamMember& team) {
        Homme::KernelVariables kv(team, num_1d_fields);
        const int ie = kv.ie;
        const int ifield = kv.iq;
        const auto ttr = Kokkos::TeamThreadRange(kv.team, 2);
        const auto tvr = Kokkos::ThreadVectorRange(kv.team, NUM_LEV);
        const int iconn_end = ucon_ptr(ie+1);
        for (int iconn = ucon_ptr(ie); iconn < iconn_end; ++iconn) {
          const auto& info = ucon(iconn);
          assert(info.kind != etoi(ConnectionSharing::MISSING));
          const int buffer_iconn = (info.sharing == etoi(ConnectionSharing::LOCAL) ?
                                    info.sharing_local_remote_iconn :
                                    iconn);
          const auto& sb = send_1d_buffers(ifield, buffer_iconn);
          const auto& f1 = fields_1d(ie, ifield);
          Kokkos::parallel_for(ttr, [&] (const int& k) {
            auto* const sbp = &sb(k, 0);
            const auto* const f1p = &f1(k, 0);
            Kokkos::parallel_for(tvr, [&] (const int& ilev) { sbp[ilev] = f1p[ilev]; });
          });
        }
      });
  }
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

  pack_min_max(m_connectivity->get_d_ucon(), m_connectivity->get_d_ucon_ptr(),
               m_1d_fields, m_send_1d_buffers, m_num_elems, m_num_1d_fields);
  Kokkos::fence();

  // ---- Send ---- //
  m_buffers_manager->sync_send_buffer(this);
  if ( ! m_send_requests.empty())
    HOMMEXX_MPI_CHECK_ERROR(MPI_Startall(m_send_requests.size(), m_send_requests.data()),
                            m_connectivity->get_comm().mpi_comm());

  // Mark send buffer as busy
  m_send_pending = true;
}

static void unpack_min_max (
  const ExecViewUnmanaged<const HaloExchangeUnstructuredConnectionInfo*> ucon,
  const ExecViewUnmanaged<const int*> ucon_ptr,
  const ExecViewUnmanaged<ExecViewManaged<Scalar[2][NUM_LEV]>**> fields_1d,
  const ExecViewUnmanaged<ExecViewUnmanaged<Scalar[2][NUM_LEV]>**> recv_1d_buffers,
  const int num_elems, const int num_1d_fields)
{
  if (OnGpu<ExecSpace>::value) {
    const ConnectionHelpers helpers;
    Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace>(0, num_elems*num_1d_fields*NUM_LEV),
      KOKKOS_LAMBDA(const int it) {
        const int ie = it / (num_1d_fields*NUM_LEV);
        const int ifield = (it / NUM_LEV) % num_1d_fields;
        const int ilev = it % NUM_LEV;
        const auto& f1 = fields_1d(ie, ifield);
        const auto iconn_beg = ucon_ptr(ie), iconn_end = ucon_ptr(ie+1);
        for (int iconn = iconn_beg; iconn < iconn_end; ++iconn) {
          const auto& r1 = recv_1d_buffers(ifield, iconn);
          f1(MIN_ID, ilev) = min(f1(MIN_ID, ilev), r1(MIN_ID, ilev));
          f1(MAX_ID, ilev) = max(f1(MAX_ID, ilev), r1(MAX_ID, ilev));
        }
      });
  } else {
    HOMMEXX_STATIC const ConnectionHelpers helpers;
    const auto num_parallel_iterations = num_elems*num_1d_fields;
    ThreadPreferences tp;
    tp.max_threads_usable = 2;
    tp.max_vectors_usable = NUM_LEV;
    const auto threads_vectors =
      DefaultThreadsDistribution<ExecSpace>::team_num_threads_vectors(
        num_parallel_iterations, tp);
    const auto policy = Kokkos::TeamPolicy<ExecSpace>(
      num_parallel_iterations, threads_vectors.first, threads_vectors.second);
    Kokkos::parallel_for(policy,
      KOKKOS_LAMBDA(const TeamMember& team) {
        Homme::KernelVariables kv(team, num_1d_fields);
        const int ie = kv.ie;
        const int ifield = kv.iq;
        const auto ttr = Kokkos::TeamThreadRange(kv.team, 2);
        const auto tvr = Kokkos::ThreadVectorRange(kv.team, NUM_LEV);
        const auto& f1 = fields_1d(ie, ifield);
        const auto iconn_beg = ucon_ptr(ie), iconn_end = ucon_ptr(ie+1);
        for (int iconn = iconn_beg; iconn < iconn_end; ++iconn) {
          const auto& r1 = recv_1d_buffers(ifield, iconn);
          Kokkos::parallel_for(ttr, [&] (const int& k) {
            const auto* const r1p = &r1(k, 0);
            auto* const f1p = &f1(k, 0);
            switch (k) {
            case MIN_ID:
              Kokkos::parallel_for(
                tvr, [&] (const int& ilev) { f1p[ilev] = min(f1p[ilev], r1p[ilev]); });
              break;
            case MAX_ID:
              Kokkos::parallel_for(
                tvr, [&] (const int& ilev) { f1p[ilev] = max(f1p[ilev], r1p[ilev]); });
              break;
            default: assert(false);
            }
          });
        }
      });
  }
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

  unpack_min_max(m_connectivity->get_d_ucon(), m_connectivity->get_d_ucon_ptr(),
                 m_1d_fields, m_recv_1d_buffers, m_num_elems, m_num_1d_fields);
  Kokkos::fence();

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

  // Check that the MpiBuffersManager is present and was setup with enough storage
  assert (m_buffers_manager);

  // Ask the buffer manager to check for reallocation and then proceed with the allocation (if needed)
  // Note: if BM already knows about our needs, and buffers were already allocated, then
  //       these two calls should not change the internal state of the BM
  m_buffers_manager->check_for_reallocation();
  m_buffers_manager->allocate_buffers();

  assert (m_3d_nlev_pack.empty() ||
          static_cast<int>(m_3d_nlev_pack.size()) == m_num_3d_fields);

  // We want to set the send/recv buffers to point to:
  //   - a portion of send/recv_buffer if info.sharing=SHARED
  //   - a portion of local_buffer if info.sharing=LOCAL
  //   - the blackhole_send/recv if info.sharing=MISSING
  // After reserving the buffer portion, update the offset by a given increment, depending on info.kind:
  //   - increment[CORNER]  = m_elem_buf_size[CORNER)] = 1  * (m_num_2d_fields + NUM_LEV*VECTOR_SIZE m_num_3d_fields)
  //   - increment[EDGE]    = m_elem_buf_size[EDGE)]   = NP * (m_num_2d_fields + NUM_LEV*VECTOR_SIZE m_num_3d_fields)
  //   - increment[MISSING] = 0 (point to the same blackhole)

  HostViewManaged<size_t[3]> h_buf_offset("");
  Kokkos::deep_copy(h_buf_offset, 0);

  // The number of Reals used in a connection on a single level:
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

  std::vector<int> slot_idx_to_elem_conn_pair, pids, pid_offsets;
  init_slot_idx_to_elem_conn_pair(slot_idx_to_elem_conn_pair, pids, pid_offsets);

  const auto& ucon = m_connectivity->get_h_ucon();
  const size_t nconn = ucon.size();
  
  m_send_1d_buffers = decltype(m_send_1d_buffers)("1d send buffer", m_num_1d_fields, nconn);
  m_recv_1d_buffers = decltype(m_recv_1d_buffers)("1d recv buffer", m_num_1d_fields, nconn);
  m_send_2d_buffers = decltype(m_send_2d_buffers)("2d send buffer", m_num_2d_fields, nconn);
  m_recv_2d_buffers = decltype(m_recv_2d_buffers)("2d recv buffer", m_num_2d_fields, nconn);
  m_send_3d_buffers = decltype(m_send_3d_buffers)("3d send buffer", m_num_3d_fields, nconn);
  m_recv_3d_buffers = decltype(m_recv_3d_buffers)("3d recv buffer", m_num_3d_fields, nconn);
  m_send_3d_int_buffers = decltype(m_send_3d_int_buffers)("3d interface send buffer", m_num_3d_int_fields, nconn);
  m_recv_3d_int_buffers = decltype(m_recv_3d_int_buffers)("3d interface recv buffer", m_num_3d_int_fields, nconn);
  const auto h_send_1d_buffers = Kokkos::create_mirror_view(m_send_1d_buffers);
  const auto h_recv_1d_buffers = Kokkos::create_mirror_view(m_recv_1d_buffers);
  const auto h_send_2d_buffers = Kokkos::create_mirror_view(m_send_2d_buffers);
  const auto h_recv_2d_buffers = Kokkos::create_mirror_view(m_recv_2d_buffers);
  const auto h_send_3d_buffers = Kokkos::create_mirror_view(m_send_3d_buffers);
  const auto h_recv_3d_buffers = Kokkos::create_mirror_view(m_recv_3d_buffers);
  const auto h_send_3d_int_buffers = Kokkos::create_mirror_view(m_send_3d_int_buffers);
  const auto h_recv_3d_int_buffers = Kokkos::create_mirror_view(m_recv_3d_int_buffers);

  ConnectionHelpers helpers;
  for (size_t k = 0; k < nconn; ++k) {
    // Map from MPI buffer index space to (elem, connection) index space.
    const auto i = slot_idx_to_elem_conn_pair[k];
    const auto& info = ucon(i);

    auto& send_buffer = h_all_send_buffers[info.sharing];
    auto& recv_buffer = h_all_recv_buffers[info.sharing];

    for (int f = 0; f < m_num_1d_fields; ++f) {
      h_send_1d_buffers(f, i) = ExecViewUnmanaged<Scalar[2][NUM_LEV]>(
        reinterpret_cast<Scalar*>(send_buffer.get() + h_buf_offset[info.sharing]));
      h_recv_1d_buffers(f, i) = ExecViewUnmanaged<Scalar[2][NUM_LEV]>(
        reinterpret_cast<Scalar*>(recv_buffer.get() + h_buf_offset[info.sharing]));
      h_buf_offset[info.sharing] += h_increment_1d[info.kind]*NUM_LEV*VECTOR_SIZE;
    }
    for (int f = 0; f < m_num_2d_fields; ++f) {
      h_send_2d_buffers(f, i) = ExecViewUnmanaged<Real*>(
        send_buffer.get() + h_buf_offset[info.sharing], helpers.CONNECTION_SIZE[info.kind]);
      h_recv_2d_buffers(f, i) = ExecViewUnmanaged<Real*>(
        recv_buffer.get() + h_buf_offset[info.sharing], helpers.CONNECTION_SIZE[info.kind]);
      h_buf_offset[info.sharing] += h_increment_2d[info.kind];
    }
    for (int f = 0; f < m_num_3d_fields; ++f) {
      const auto nlev_3d = m_3d_nlev_pack.empty() ? NUM_LEV : m_3d_nlev_pack[f];
      h_send_3d_buffers(f, i) = ExecViewUnmanaged<Scalar**>(
        reinterpret_cast<Scalar*>(send_buffer.get() + h_buf_offset[info.sharing]),
        helpers.CONNECTION_SIZE[info.kind], nlev_3d);
      h_recv_3d_buffers(f, i) = ExecViewUnmanaged<Scalar**>(
        reinterpret_cast<Scalar*>(recv_buffer.get() + h_buf_offset[info.sharing]),
        helpers.CONNECTION_SIZE[info.kind], nlev_3d);
      h_buf_offset[info.sharing] += h_increment_3d[info.kind]*nlev_3d*VECTOR_SIZE;
    }
    for (int f = 0; f < m_num_3d_int_fields; ++f) {
      h_send_3d_int_buffers(f, i) = ExecViewUnmanaged<Scalar**>(
        reinterpret_cast<Scalar*>(send_buffer.get() + h_buf_offset[info.sharing]),
        helpers.CONNECTION_SIZE[info.kind], NUM_LEV_P);
      h_recv_3d_int_buffers(f, i) = ExecViewUnmanaged<Scalar**>(
        reinterpret_cast<Scalar*>(recv_buffer.get() + h_buf_offset[info.sharing]),
        helpers.CONNECTION_SIZE[info.kind], NUM_LEV_P);
      h_buf_offset[info.sharing] += h_increment_3d[info.kind]*NUM_LEV_P*VECTOR_SIZE;
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
  // Sanity check: compute the buffers sizes for this boundary exchange, and
  // check that the final offsets match them.
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
    const auto mpi_comm = m_connectivity->get_comm().mpi_comm();
    const size_t npids = pids.size();
    free_requests();
    m_send_requests.resize(npids);
    m_recv_requests.resize(npids);
    MPIViewManaged<Real*>::pointer_type send_ptr = buffers_manager->get_mpi_send_buffer().data();
    MPIViewManaged<Real*>::pointer_type recv_ptr = buffers_manager->get_mpi_recv_buffer().data();
    int offset = 0;
    for (size_t ip = 0; ip < npids; ++ip) {
      int count = 0;
      for (int k = pid_offsets[ip]; k < pid_offsets[ip+1]; ++k) {
        const auto i = slot_idx_to_elem_conn_pair[k];
        const auto& info = ucon(i);
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
    size_t i, ord;
    int pid;
    bool operator< (const IP& o) const {
      if (pid < o.pid) return true;
      if (pid > o.pid) return false;
      return ord < o.ord;
    }
  };

  const auto& ucon = m_connectivity->get_h_ucon();
  const auto& ucon_ptr = m_connectivity->get_h_ucon_ptr();
  const size_t nconn = ucon.size();
  const int mce = m_connectivity->get_max_corner_elements();
  const int n_idx_per_elem = 8*mce;
  std::vector<IP> i2remote(nconn);

  for (size_t k = 0; k < nconn; ++k) {
    const auto& info = ucon(k);
    auto& i2r = i2remote[k];
    // Original sequence through connections.
    i2r.i = k;
    // An ordering of the message buffer upon which both members of the
    // communication pair agree.
    const auto& lgp = info.local.gid < info.remote.gid ? info.local : info.remote;
    i2r.ord = lgp.gid*n_idx_per_elem + lgp.dir*mce + lgp.dir_idx;
    // If local, indicate with -1, which is < the smallest pid of 0.
    i2r.pid = -1;
    if (info.sharing == etoi(ConnectionSharing::SHARED))
      i2r.pid = info.remote_pid;
  }

  // Sort so that, first, all connections having the same remote_pid are
  // contiguous; second, within such a block, ord is ascending. The first lets
  // us set up comm buffers so monolithic messages slot right in. The second
  // means that the send and recv partners agree on how a monolithic message is
  // packed.
  std::sort(i2remote.begin(), i2remote.end());

  // Collect the unique remote_pids and get the offsets of the contiguous blocks
  // of them.
  slot_idx_to_elem_conn_pair.resize(nconn);
  pids.clear();
  pid_offsets.clear();
  int prev_pid = -2;
  for (size_t k = 0; k < nconn; ++k) {
    const auto& i2r = i2remote[k];
    if (i2r.pid > prev_pid && i2r.pid != -1) {
      pids.push_back(i2r.pid);
      pid_offsets.push_back(k);
      prev_pid = i2r.pid;
    }
    slot_idx_to_elem_conn_pair[k] = i2r.i;
  }
  pid_offsets.push_back(nconn);
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
