#include "compose_slmm_islmpi.hpp"
#include "compose_slmm_islmpi_buf.hpp"

namespace homme {
namespace islmpi {

#ifdef COMPOSE_PORT
/* GPU metadata are arranged differently than described below. The scheme is the
   following:
        (x-bulk-data-offset int
         #x-in-rank         i
         (lid-on-rank       i     only packed if #x in lid > 0
          lev               short
          #x)               s
             *#x-in-rank)
*/

struct Accum {
  Int mos, sendcount, xos, qos;
  SLMM_KIF Accum () : mos(0), sendcount(0), xos(0), qos(0) {}
  SLMM_KIF void operator+= (const volatile Accum& o) volatile {
    mos += o.mos; sendcount += o.sendcount; xos += o.xos; qos += o.qos;
  }
};

template <typename MT>
void pack_dep_points_sendbuf_pass1_scan (IslMpi<MT>& cm, const bool trajectory) {
  ko::fence();
  deep_copy(cm.nx_in_rank_h, cm.nx_in_rank);
  const auto& sendbufs = cm.sendbuf;
  const auto& lid_on_ranks = cm.lid_on_rank;
  const auto& nx_in_rank = cm.nx_in_rank;
  const auto& nx_in_lids = cm.nx_in_lid;
  const auto& x_bulkdata_offset = cm.x_bulkdata_offset;
  const auto& sendcounts = cm.sendcount;
  const auto& blas = cm.bla;
  const auto nlev = cm.nlev;
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  const Int ndim = trajectory ? cm.dep_points_ndim : 3;
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const Int lid_on_rank_n = cm.lid_on_rank_h(ri).n();
    const auto f = COMPOSE_LAMBDA (const int idx, Accum& a, const bool fin) {
      auto&& sendbuf = sendbufs(ri);
      if (idx == 0) {
        const auto cnt = setbuf(sendbuf, 0, 0, 0, fin);
        a.mos += cnt;
        a.sendcount += cnt;
      }
      if (nx_in_rank(ri) == 0) return;
      const auto&& lid_on_rank = lid_on_ranks(ri);
      const Int lidi = idx / nlev;
      const Int lev = idx % nlev;
      const auto nx_in_lid = nx_in_lids(ri,lidi);
      if (nx_in_lid == 0) return;
      auto&& bla = blas(ri);
      auto& t = bla(lidi,lev);
      slmm_kernel_assert_high(t.cnt == 0);
      const Int nx = t.xptr;
      if (fin) {
        t.qptr = a.qos;
        if (nx == 0) {
          t.xptr = -1;
          return;
        }
      }
      if (nx > 0) {
        const auto dos = setbuf(sendbuf, a.mos, lid_on_rank(lidi), lev, nx, fin);
        a.mos += dos;
        a.sendcount += dos + ndim*nx;
        if (fin) t.xptr = a.xos;
        a.xos += ndim*nx;
        a.qos += trajectory ? nx : 2 + nx;
      }
    };
    Accum a;
    ko::parallel_scan(ko::RangePolicy<typename MT::DES>(0, lid_on_rank_n*nlev), f, a);
    const auto g = COMPOSE_LAMBDA (const int) {
      auto&& sendbuf = sendbufs(ri);
      setbuf(sendbuf, 0, a.mos /* offset to x bulk data */, nx_in_rank(ri));
      x_bulkdata_offset(ri) = a.mos;
      sendcounts(ri) = a.sendcount;
    };
    ko::parallel_for(ko::RangePolicy<typename MT::DES>(0, 1), g);
  }
  ko::fence();
  deep_copy(cm.sendcount_h, cm.sendcount);
}
#endif

/* Pack the departure points (x). We use two passes. We also set up the q
   metadata. Two passes let us do some efficient tricks that are not available
   with one pass. Departure point and q messages are formatted as follows:
    xs: (x-bulk-data-offset int                                     <-
         #x-in-rank         i                                        |
         (lid-on-rank       i     only packed if #x in lid > 0       |
          #x-in-lid         i     > 0                                |- meta data
          (lev              i     only packed if #x in (lid,lev) > 0 |
           #x)              i     > 0                                |
              *#lev) *#lid                                          <-
         x                  3 real                                  <-- bulk data
          *#x-in-rank) *#rank
    qs: (q-extrema    2 qsize r   (min, max) packed together
         q              qsize r
          *#x) *#lev *#lid *#rank
 */
template <typename MT>
void pack_dep_points_sendbuf_pass1_noscan (IslMpi<MT>& cm, const bool trajectory) {
#ifdef COMPOSE_PORT
  ko::fence();
  deep_copy(cm.nx_in_rank_h, cm.nx_in_rank);
  deep_copy(cm.nx_in_lid_h, cm.nx_in_lid);
  deep_copy(cm.bla_h, cm.bla);
#endif
  const Int nrmtrank = static_cast<Int>(cm.ranks.size()) - 1;
  const Int ndim = trajectory ? cm.dep_points_ndim : 3;
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp for
#endif
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    auto&& sendbuf = cm.sendbuf_meta_h(ri);
    const auto&& lid_on_rank = cm.lid_on_rank_h(ri);
    // metadata offset, x bulk data offset, q bulk data offset
    Int mos = 0, xos = 0, qos = 0, sendcount = 0, cnt;
    cnt = setbuf(sendbuf, mos, 0, 0); // empty space for later
    mos += cnt;
    sendcount += cnt;
    if (cm.nx_in_rank_h(ri) == 0) {
      setbuf(sendbuf, 0, mos, 0);
      cm.x_bulkdata_offset_h(ri) = mos;
      cm.sendcount_h(ri) = sendcount;
      continue;
    }
    auto&& bla = cm.bla_h(ri);
    for (Int lidi = 0, lidn = lid_on_rank.n(); lidi < lidn; ++lidi) {
      auto nx_in_lid = cm.nx_in_lid_h(ri,lidi);
      if (nx_in_lid == 0) continue;
      cnt = setbuf(sendbuf, mos, lid_on_rank(lidi), nx_in_lid);
      mos += cnt;
      sendcount += cnt;
      for (Int lev = 0; lev < cm.nlev; ++lev) {
        auto& t = bla(lidi,lev);
        t.qptr = qos;
        slmm_assert_high(t.cnt == 0);
        const Int nx = t.xptr;
        if (nx == 0) {
          t.xptr = -1;
          continue;
        }
        slmm_assert_high(nx > 0);
        const auto dos = setbuf(sendbuf, mos, lev, nx);
        mos += dos;
        sendcount += dos + ndim*nx;
        t.xptr = xos;
        xos += ndim*nx;
        qos += trajectory ? nx : 2 + nx;
        nx_in_lid -= nx;
      }
      slmm_assert(nx_in_lid == 0);
    }
    setbuf(sendbuf, 0, mos /* offset to x bulk data */, cm.nx_in_rank_h(ri));
    cm.x_bulkdata_offset_h(ri) = mos;
    cm.sendcount_h(ri) = sendcount;
  }
#ifdef COMPOSE_PORT
  deep_copy(cm.sendcount, cm.sendcount_h);
  deep_copy(cm.x_bulkdata_offset, cm.x_bulkdata_offset_h);
  deep_copy(cm.bla, cm.bla_h);  
#endif
#ifdef COMPOSE_PORT_SEPARATE_VIEWS
  // Copy metadata chunks to device sendbuf.
  slmm_assert(cm.sendbuf.n() == nrmtrank);
  for (Int ri = 0; ri < nrmtrank; ++ri) {
    const auto n = cm.x_bulkdata_offset_h(ri);
    assert(n <= cm.sendmetasz[ri]);
    if (n > 0)
      ko::deep_copy(ko::View<Real*, typename MT::DES>(cm.sendbuf.get_h(ri).data(), n),
                    ko::View<Real*, typename MT::HES>(cm.sendbuf_meta_h(ri).data(), n));
  }
#endif
}

template <typename MT>
void pack_dep_points_sendbuf_pass1 (IslMpi<MT>& cm, const bool trajectory) {
#if defined COMPOSE_PORT && ! defined COMPOSE_PACK_NOSCAN
  if (ko::OnGpu<typename MT::DES>::value)
    pack_dep_points_sendbuf_pass1_scan(cm, trajectory);
  else
#endif
    pack_dep_points_sendbuf_pass1_noscan(cm, trajectory);
}

template <typename MT>
void pack_dep_points_sendbuf_pass2 (IslMpi<MT>& cm, const DepPoints<MT>& dep_points,
                                    const bool trajectory) {
  const auto myrank = cm.p->rank();
#ifdef COMPOSE_PORT
  const Int start = 0, end = cm.mylid_with_comm_h.n();
#else
  const int tid = get_tid();
  ConstExceptGnu Int
    start = cm.mylid_with_comm_tid_ptr_h(tid),
    end = cm.mylid_with_comm_tid_ptr_h(tid+1);
#endif
  {
    auto ed = cm.ed_d.unmanaged();
    const auto mylid_with_comm = cm.mylid_with_comm_d.unmanaged();
    ko::parallel_for(
      ko::RangePolicy<typename MT::DES>(start, end),
      COMPOSE_LAMBDA (const Int& ptr) {
        const Int tci = mylid_with_comm(ptr);
        ed(tci).rmt.clear();
      });
  }
  {
    ConstExceptGnu Int np2 = cm.np2, nlev = cm.nlev, qsize = cm.qsize;
    ConstExceptGnu Int ndim = trajectory ? cm.dep_points_ndim : 3;
    const auto& ed_d = cm.ed_d;
    const auto& mylid_with_comm_d = cm.mylid_with_comm_d;
    const auto& sendbuf = cm.sendbuf;
    const auto& x_bulkdata_offset = cm.x_bulkdata_offset;
    const auto& bla = cm.bla;
#ifdef COMPOSE_HORIZ_OPENMP
    const auto horiz_openmp = cm.horiz_openmp;
    const auto& ri_lidi_locks = cm.ri_lidi_locks;
#endif
    const auto f = COMPOSE_LAMBDA (const Int& ki) {
      const Int ptr = start + ki/(nlev*np2);
      const Int lev = (ki/np2) % nlev;
      const Int k = ki % np2;
      const Int tci = mylid_with_comm_d(ptr);
      auto& ed = ed_d(tci);
      const Int sci = ed.src(lev,k);
      const auto& nbr = ed.nbrs(sci);
      if (nbr.rank == myrank) return;
      const Int ri = nbr.rank_idx;
      const Int lidi = nbr.lid_on_rank_idx;
      auto&& sb = sendbuf(ri);
#ifdef COMPOSE_HORIZ_OPENMP
      omp_lock_t* lock;
      if (horiz_openmp) {
        lock = &ri_lidi_locks(ri,lidi);
        omp_set_lock(lock);
      }
#endif
      Int xptr, qptr, cnt; {
        auto& t = bla(ri,lidi,lev);
#ifdef COMPOSE_PORT
        cnt = ko::atomic_fetch_add(static_cast<volatile Int*>(&t.cnt), 1);
#else
        cnt = t.cnt;
        ++t.cnt;
#endif
        qptr = t.qptr;
        xptr = x_bulkdata_offset(ri) + t.xptr + ndim*cnt;
      }
#ifdef COMPOSE_HORIZ_OPENMP
      if (horiz_openmp) omp_unset_lock(lock);
#endif
      slmm_kernel_assert_high(xptr > 0);
      for (Int i = 0; i < ndim; ++i)
        sb(xptr + i) = dep_points(tci,lev,k,i);
      auto& item = ed.rmt.atomic_inc_and_return_next();
      if (trajectory) {
        item.q_extrema_ptr = item.q_ptr = ndim*(qptr + cnt);
      } else {
        item.q_extrema_ptr = qsize * qptr;
        item.q_ptr = item.q_extrema_ptr + qsize*(2 + cnt);
      }
      item.lev = lev;
      item.k = k;
    };
    ko::fence();
    ko::parallel_for(
      ko::RangePolicy<typename MT::DES>(0, (end - start)*nlev*np2), f);
    ko::fence();
  }
}

template void pack_dep_points_sendbuf_pass1(
  IslMpi<ko::MachineTraits>& cm, const bool trajectory);
template void pack_dep_points_sendbuf_pass2(
  IslMpi<ko::MachineTraits>& cm, const DepPoints<ko::MachineTraits>& dep_points,
  const bool trajectory);

} // namespace islmpi
} // namespace homme
