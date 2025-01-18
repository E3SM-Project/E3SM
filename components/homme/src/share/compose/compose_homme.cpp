#include "compose_homme.hpp"

namespace homme {

template <typename MT>
TracerArrays<MT>::TracerArrays (Int nelemd_, Int nlev_, Int np_, Int qsize_, Int qsized_)
  : nelemd(nelemd_), nlev(nlev_), np(np_), np2(np*np), qsize(qsize_), qsized(qsized_),
    pspheremp(nelemd, np2),
    pdp(nelemd, np2, nlev), pdp3d(nelemd, np2, nlev, -1, 3),
    pqdp(nelemd, np2, nlev, qsized, 2), pq(nelemd, np2, nlev, qsized),
#if defined COMPOSE_PORT
    q_min("q_min", nelemd, qsize, np2, nlev),
    q_max("q_max", nelemd, qsize, np2, nlev)
#else
    spheremp(pspheremp), dp(pdp), dp3d(pdp3d), qdp(pqdp), q(pq)
#endif
{}

#if defined COMPOSE_PORT
template <typename MT>
void TracerArrays<MT>::alloc_if_not () {
  if (spheremp.size() > 0) return;
  spheremp = decltype(spheremp)("spheremp", nelemd, np2);
  dp = decltype(dp)("dp", nelemd, np2, nlev);
  dp3d = decltype(dp3d)("dp3d", nelemd, 3, np2, nlev);
  qdp = decltype(qdp)("qdp", nelemd, 2, qsize, np2, nlev);
  q = decltype(q)("q", nelemd, qsize, np2, nlev);
  assert(dep_points.size() > 0);
}
#endif

template <typename MT>
void sl_traj_h2d (TracerArrays<MT>& ta, Real* dep_points, Real* vnode,
                  Real* vdep, Int ndim) {
#if defined COMPOSE_PORT
# if defined COMPOSE_HORIZ_OPENMP
# pragma omp master
  {
# endif
  ko::fence();
  ta.alloc_if_not();
  const Int nelemd = ta.nelemd, qsize = ta.qsize, np2 = ta.np2, nlev = ta.nlev;
  const DepPointsH<MT> cart_h(dep_points, nelemd, nlev, np2, ndim);
  ko::deep_copy(ta.dep_points, cart_h);
  if (vnode) {
    const DepPointsH<MT> h(vnode, nelemd, nlev, np2, ndim);
    ko::deep_copy(ta.vnode, h);
  }
  if (vdep) {
    const DepPointsH<MT> h(vdep, nelemd, nlev, np2, ndim);
    ko::deep_copy(ta.vdep, h);    
  }
# ifdef COMPOSE_HORIZ_OPENMP
  }
# pragma omp barrier
# endif
#endif
}

template <typename MT>
void sl_traj_d2h (const TracerArrays<MT>& ta, Real* dep_points, Real* vnode,
                  Real* vdep, Int ndim) {
#if defined COMPOSE_PORT
# if defined COMPOSE_HORIZ_OPENMP
# pragma omp master
  {
# endif
  ko::fence();
  const auto q_m = ko::create_mirror_view(ta.q);
  const Int nelemd = ta.nelemd, np2 = ta.np2, nlev = ta.nlev;
  const DepPointsH<MT> dep_points_h(dep_points, nelemd, nlev, np2, ndim);
  ko::deep_copy(dep_points_h, ta.dep_points);
  if (vnode) {
    const DepPointsH<MT> h(vnode, nelemd, nlev, np2, ndim);
    ko::deep_copy(h, ta.vnode);
  }
  if (vdep) {
    const DepPointsH<MT> h(vdep, nelemd, nlev, np2, ndim);
    ko::deep_copy(h, ta.vdep);    
  }
# ifdef COMPOSE_HORIZ_OPENMP
  }
# pragma omp barrier
# endif
#endif  
}

template <typename MT>
void sl_h2d (TracerArrays<MT>& ta, bool transfer, Real* dep_points, Int ndim) {
#if defined COMPOSE_PORT
# if defined COMPOSE_HORIZ_OPENMP
# pragma omp master
  {
# endif
  ko::fence();
  ta.alloc_if_not();
  const Int nelemd = ta.nelemd, qsize = ta.qsize, np2 = ta.np2, nlev = ta.nlev;
  const DepPointsH<MT> cart_h(dep_points, nelemd, nlev, np2, ndim);
  ko::deep_copy(ta.dep_points, cart_h);
  if (transfer) {
    const auto qdp_m = ko::create_mirror_view(ta.qdp);
    const auto dp_m = ko::create_mirror_view(ta.dp);
    const auto q_m = ko::create_mirror_view(ta.q);
    for (Int ie = 0; ie < nelemd; ++ie)
      for (Int iq = 0; iq < qsize; ++iq)
        for (Int k = 0; k < np2; ++k)
          for (Int lev = 0; lev < nlev; ++lev) {
            for (Int qtl = 0; qtl < 2; ++qtl)
              qdp_m(ie,qtl,iq,k,lev) = ta.pqdp(ie,qtl,iq,k,lev);
            q_m(ie,iq,k,lev) = ta.pq(ie,iq,k,lev);
          }
    for (Int ie = 0; ie < nelemd; ++ie)
      for (Int k = 0; k < np2; ++k)
        for (Int lev = 0; lev < nlev; ++lev)
          dp_m(ie,k,lev) = ta.pdp(ie,k,lev);
    ko::deep_copy(ta.qdp, qdp_m);
    ko::deep_copy(ta.dp, dp_m);
    ko::deep_copy(ta.q, q_m);
  }
# ifdef COMPOSE_HORIZ_OPENMP
  }
# pragma omp barrier
# endif
#endif
}

template <typename MT>
void sl_d2h (const TracerArrays<MT>& ta, bool transfer, Real* dep_points, Int ndim,
             Real* minq, Real* maxq) {
#if defined COMPOSE_PORT
  if ( ! transfer) return;
# if defined COMPOSE_HORIZ_OPENMP
# pragma omp master
  {
# endif
  ko::fence();
  const auto q_m = ko::create_mirror_view(ta.q);
  const Int nelemd = ta.nelemd, qsize = ta.qsize, np2 = ta.np2, nlev = ta.nlev;
  ko::deep_copy(q_m, ta.q);
  for (Int ie = 0; ie < nelemd; ++ie)
    for (Int iq = 0; iq < qsize; ++iq)
      for (Int k = 0; k < np2; ++k)
        for (Int lev = 0; lev < nlev; ++lev)
          ta.pq(ie,iq,k,lev) = q_m(ie,iq,k,lev);
  const DepPointsH<MT> dep_points_h(dep_points, nelemd, nlev, np2, ndim);
  const QExtremaH<MT>
    q_min_h(minq, nelemd, qsize, np2, nlev),
    q_max_h(maxq, nelemd, qsize, np2, nlev);
  ko::deep_copy(dep_points_h, ta.dep_points);
  ko::deep_copy(q_min_h, ta.q_min);
  ko::deep_copy(q_max_h, ta.q_max);
# ifdef COMPOSE_HORIZ_OPENMP
  }
# pragma omp barrier
# endif
#endif  
}

template <typename MT>
void cedr_h2d (const TracerArrays<MT>& ta, bool transfer) {
#if defined COMPOSE_PORT
  if ( ! transfer) return;
# if defined COMPOSE_HORIZ_OPENMP
# pragma omp master
  {
# endif
  ko::fence();
  const auto dp3d_m = ko::create_mirror_view(ta.dp3d);
  const auto q_m = ko::create_mirror_view(ta.q);
  const auto spheremp_m = ko::create_mirror_view(ta.spheremp);
  const Int nelemd = ta.nelemd, qsize = ta.qsize, np2 = ta.np2, nlev = ta.nlev, np1 = ta.np1;
  for (Int ie = 0; ie < nelemd; ++ie)
    for (Int iq = 0; iq < qsize; ++iq)
      for (Int k = 0; k < np2; ++k)
        for (Int lev = 0; lev < nlev; ++lev)
          q_m(ie,iq,k,lev) = ta.pq(ie,iq,k,lev);
  for (Int ie = 0; ie < nelemd; ++ie)
    for (Int k = 0; k < np2; ++k)
      for (Int lev = 0; lev < nlev; ++lev)
        dp3d_m(ie,np1,k,lev) = ta.pdp3d(ie,np1,k,lev);
  for (Int ie = 0; ie < nelemd; ++ie)
    for (Int k = 0; k < np2; ++k)
      spheremp_m(ie,k) = ta.pspheremp(ie,k);
  ko::deep_copy(ta.dp3d, dp3d_m);
  ko::deep_copy(ta.q, q_m);
  ko::deep_copy(ta.spheremp, spheremp_m);
# ifdef COMPOSE_HORIZ_OPENMP
  }
# pragma omp barrier
# endif
#endif  
}

template <typename MT>
void cedr_d2h (const TracerArrays<MT>& ta, bool transfer) {
#if defined COMPOSE_PORT
  if ( ! transfer) return;
# if defined COMPOSE_HORIZ_OPENMP
# pragma omp master
  {
# endif
  ko::fence();
  const auto q_m = ko::create_mirror_view(ta.q);
  const auto qdp_m = ko::create_mirror_view(ta.qdp);
  const Int nelemd = ta.nelemd, qsize = ta.qsize, np2 = ta.np2, nlev = ta.nlev,
    n1_qdp = ta.n1_qdp;
  ko::deep_copy(qdp_m, ta.qdp);
  ko::deep_copy(q_m, ta.q);
  for (Int ie = 0; ie < nelemd; ++ie)
    for (Int iq = 0; iq < qsize; ++iq)
      for (Int k = 0; k < np2; ++k)
        for (Int lev = 0; lev < nlev; ++lev) {
          ta.pqdp(ie,n1_qdp,iq,k,lev) = qdp_m(ie,n1_qdp,iq,k,lev);
          ta.pq(ie,iq,k,lev) = q_m(ie,iq,k,lev);
        }
# ifdef COMPOSE_HORIZ_OPENMP
  }
# pragma omp barrier
# endif
#endif
}

TracerArrays<ko::MachineTraits>::Ptr& get_instance () {
  static typename TracerArrays<ko::MachineTraits>::Ptr p;
  return p;
}

TracerArrays<ko::MachineTraits>::Ptr
init_tracer_arrays (Int nelemd, Int nlev, Int np2, Int qsize, Int qsized) {
  auto& p = get_instance();
  if (p == nullptr)
    p = std::make_shared<TracerArrays<ko::MachineTraits> >(nelemd, nlev, np2, qsize, qsized);
  return p;
}

TracerArrays<ko::MachineTraits>::Ptr get_tracer_arrays () {
  return get_instance();
}

void delete_tracer_arrays () {
  auto& p = get_instance();
  p = nullptr;
}

template struct TracerArrays<ko::MachineTraits>;
template void sl_traj_h2d(TracerArrays<ko::MachineTraits>& ta,
                          Real*, Real*, Real*, Int ndim);
template void sl_traj_d2h(const TracerArrays<ko::MachineTraits>& ta,
                          Real*, Real*, Real*, Int ndim);
template void sl_h2d(TracerArrays<ko::MachineTraits>& ta, bool transfer,
                     Real* dep_points, Int ndim);
template void sl_d2h(const TracerArrays<ko::MachineTraits>& ta, bool transfer,
                     Real* dep_points, Int ndim, Real* minq, Real* maxq);
template void cedr_h2d(const TracerArrays<ko::MachineTraits>& ta, bool transfer);
template void cedr_d2h(const TracerArrays<ko::MachineTraits>& ta, bool transfer);

} // namespace homme
