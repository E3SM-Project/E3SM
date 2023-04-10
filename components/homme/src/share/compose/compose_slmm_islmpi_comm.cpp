#include "compose_slmm_islmpi.hpp"

namespace homme {
namespace islmpi {
// mylid_with_comm(rankidx) is a list of element LIDs that have relations with
// other elements on other ranks. For horizontal threading, need to find the
// subsets that fit within the usual horizontal-threading nets:nete ranges.
template <typename MT>
void init_mylid_with_comm_threaded (IslMpi<MT>& cm, const Int& nets, const Int& nete) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    const int nthr = get_num_threads();
#ifndef COMPOSE_PORT
    cm.rwork = typename IslMpi<MT>::template ArrayD<Real**>("rwork", nthr, cm.qsize);
#endif
    cm.mylid_with_comm_tid_ptr_h.reset_capacity(nthr+1, true);
    cm.horiz_openmp = get_num_threads() > 1;
  }
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
  const int tid = get_tid();
  const auto& beg = std::lower_bound(cm.mylid_with_comm_h.begin(),
                                     cm.mylid_with_comm_h.end(), nets);
  slmm_assert(cm.p->size() == 1 || beg != cm.mylid_with_comm_h.end());
  cm.mylid_with_comm_tid_ptr_h(tid) = beg - cm.mylid_with_comm_h.begin();
  if (tid == cm.mylid_with_comm_tid_ptr_h.n() - 2) {
    const auto& end = std::lower_bound(cm.mylid_with_comm_h.begin(),
                                       cm.mylid_with_comm_h.end(), nete+1);
    cm.mylid_with_comm_tid_ptr_h(tid+1) = end - cm.mylid_with_comm_h.begin();
  }
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

template <typename MT>
void setup_irecv (IslMpi<MT>& cm, const bool skip_if_empty) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp master
#endif
  {
    const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
    cm.recvreq.clear();
    for (Int ri = 0, nri = 0; ri < nrmtrank; ++ri) {
      if (skip_if_empty && cm.nx_in_rank_h(ri) == 0) continue;
      // The count is just the number of slots available, which can be larger
      // than what is actually being received.
      cm.recvreq_ri(nri++) = ri;
      cm.recvreq.inc();
#ifdef COMPOSE_MPI_ON_HOST
      auto&& recvbuf = cm.recvbuf_h(ri);
#else
      auto&& recvbuf = cm.recvbuf.get_h(ri);
#endif
      mpi::irecv(*cm.p, recvbuf.data(), recvbuf.n(), cm.ranks(ri), 42,
                 &cm.recvreq.back());
    }
  }
}

template <typename MT>
void isend (IslMpi<MT>& cm, const bool want_req, const bool skip_if_empty) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
    for (Int ri = 0; ri < nrmtrank; ++ri) {
      if (skip_if_empty && cm.sendcount_h(ri) == 0) continue;
#ifdef COMPOSE_MPI_ON_HOST
      auto&& sendbuf = cm.sendbuf_h(ri);
      typedef typename IslMpi<MT>::template ArrayH<Real*> ArrayH;
      typedef typename IslMpi<MT>::template ArrayD<Real*> ArrayD;
      Kokkos::deep_copy(ArrayH(sendbuf.data(), cm.sendcount_h(ri)),
                        ArrayD(cm.sendbuf.get_h(ri).data(), cm.sendcount_h(ri)));
#else
      auto&& sendbuf = cm.sendbuf.get_h(ri);
#endif
      mpi::isend(*cm.p, sendbuf.data(), cm.sendcount_h(ri),
                 cm.ranks(ri), 42, want_req ? &cm.sendreq(ri) : nullptr);
    }
  }
}

template <typename MT>
void wait_on_send (IslMpi<MT>& cm, const bool skip_if_empty) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp master
#endif
  {
    for (Int ri = 0; ri < cm.sendreq.n(); ++ri) {
      if (skip_if_empty && cm.sendcount_h(ri) == 0) continue;
      mpi::wait(&cm.sendreq(ri));
    }
  }
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

template <typename MT>
void wait_on_recv (IslMpi<MT>& cm) {
#ifdef COMPOSE_MPI_ON_HOST
  typedef typename IslMpi<MT>::template ArrayH<Real*> ArrayH;
  typedef typename IslMpi<MT>::template ArrayD<Real*> ArrayD;
  const int nreq = cm.recvreq.n();
  for (Int i = 0; i < nreq; ++i) {
    Int reqi;
    MPI_Status stat;
    mpi::waitany(nreq, cm.recvreq.data(), &reqi, &stat);
    const Int ri = cm.recvreq_ri(reqi);
    int count;
    MPI_Get_count(&stat, mpi::get_type<Real>(), &count);
    Kokkos::deep_copy(ArrayD(cm.recvbuf.get_h(ri).data(), count),
                      ArrayH(cm.recvbuf_h(ri).data(), count));
  }
#else
  mpi::waitall(cm.recvreq.n(), cm.recvreq.data());
#endif
}

template <typename MT>
void recv_and_wait_on_send (IslMpi<MT>& cm) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp master
#endif
  {
    mpi::waitall(cm.sendreq.n(), cm.sendreq.data());
    wait_on_recv(cm);
  }
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

template <typename MT>
void recv (IslMpi<MT>& cm, const bool skip_if_empty) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp master
#endif
  {
    wait_on_recv(cm);
  }
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

template void init_mylid_with_comm_threaded(
  IslMpi<ko::MachineTraits>& cm, const Int& nets, const Int& nete);
template void setup_irecv(IslMpi<ko::MachineTraits>& cm, const bool skip_if_empty);
template void isend(IslMpi<ko::MachineTraits>& cm, const bool want_req,
                    const bool skip_if_empty);
template void recv_and_wait_on_send(IslMpi<ko::MachineTraits>& cm);
template void wait_on_send(IslMpi<ko::MachineTraits>& cm, const bool skip_if_empty);
template void recv(IslMpi<ko::MachineTraits>& cm, const bool skip_if_empty);

} // namespace islmpi
} // namespace homme
