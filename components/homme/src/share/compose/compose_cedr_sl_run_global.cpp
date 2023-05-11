#include "compose_cedr.hpp"
#include "compose_cedr_cdr.hpp"
#include "compose_cedr_sl.hpp"

namespace homme {
namespace sl {

template <typename MT, typename Ie2gci, typename Qdp, typename Dp3d>
ko::EnableIfNotOnGpu<MT> warn_on_Qm_prev_negative (
  Real Qm_prev, Int rank, Int ie, const Ie2gci& ie2gci, Int np2, Int spli,
  Int k0, Int q, Int ti, Int sbli, Int lci, Int k, Int n0_qdp, Int np1,
  const Qdp& qdp_p, const Dp3d& dp3d_c)
{
  static bool first = true;
  if (first) {
    first = false;
    std::stringstream ss;
    ss << "Qm_prev < -0.5: Qm_prev = " << Qm_prev
       << " on rank " << rank
       << " at (ie,gid,spli,k0,q,ti,sbli,lci,k,n0_qdp,tl_np1) = ("
       << ie << "," << ie2gci[ie] << "," << spli << "," << k0 << ","
       << q << "," << ti << "," << sbli << "," << lci << "," << k << ","
       << n0_qdp << "," << np1 << ")\n";
    ss << "Qdp(:,:,k,q,n0_qdp) = [";
    for (Int g = 0; g < np2; ++g)
      ss << " " << qdp_p(ie,n0_qdp,q,g,k);
    ss << "]\n";
    ss << "dp3d(:,:,k,tl_np1) = [";
    for (Int g = 0; g < np2; ++g)
      ss << " " << dp3d_c(ie,np1,g,k);
    ss << "]\n";
    pr(ss.str());
  }
}

template <typename MT, typename Ie2gci, typename Qdp, typename Dp3d>
COMPOSE_FUNCTION ko::EnableIfOnGpu<MT> warn_on_Qm_prev_negative (
  Real Qm_prev, Int rank, Int ie, const Ie2gci& ie2gci, Int np2, Int spli,
  Int k0, Int q, Int ti, Int sbli, Int lci, Int k, Int n0_qdp, Int np1,
  const Qdp& qdp_p, const Dp3d& dp3d_c)
{}

template <typename MT>
static void run_cdr (CDR<MT>& q) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
  q.cdr->run();
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

template <int np_, typename MT, typename CDRT>
void run_global (CDR<MT>& cdr, CDRT* cedr_cdr_p,
                 const Data& d, Real* q_min_r, const Real* q_max_r,
                 const Int nets, const Int nete) {
  Timer t("01_write_global");
  const auto& ta = *d.ta;
  cedr_assert(ta.np == np_);
  static const Int np2 = np_*np_;
  const Int nlev = ta.nlev, qsize = ta.qsize, nlevwrem = cdr.nsuplev*cdr.nsublev;
  cedr_assert(np_ <= 4);
  
#ifdef COMPOSE_PORT
  const auto& q_min = ta.q_min;
  const auto& q_max = ta.q_max;
#else
  const QExtremaH<MT>      q_min(q_min_r, ta.nelemd, ta.qsize, ta.nlev, ta.np2);
  const QExtremaHConst<MT> q_max(q_max_r, ta.nelemd, ta.qsize, ta.nlev, ta.np2);
#endif
  const auto np1 = ta.np1;
  const auto n0_qdp = ta.n0_qdp;
  const auto& spheremp = ta.spheremp;
  const auto& dp3d_c = ta.dp3d;
  const auto& qdp_p = ta.qdp;
  const auto& q_c = ta.q;

  const Int nsublev = cdr.nsublev;
  const Int nsuplev = cdr.nsuplev;
  const auto rank = cdr.p->rank();
  const auto cdr_over_super_levels = cdr.cdr_over_super_levels;
  const auto caas_in_suplev = cdr.caas_in_suplev;
  const auto& nonnegs = cdr.nonneg;
  const auto& ie2lci = cdr.ie2lci;
  const auto& ie2gci = cdr.ie2gci;
  const typename CDRT::DeviceOp
#ifndef COMPOSE_PORT
    & // When running F90 Homme's threading scheme, we can't create new views.
#endif
    cedr_cdr = cedr_cdr_p->get_device_op();
  if (cedr::impl::OnGpu<typename MT::DES>::value) {
    const Int n = ta.nelemd*nlev*qsize*np2;
    ko::View<Real*> q_min_1d(q_min.data(), n);
    ko::parallel_for(ko::RangePolicy<typename MT::DES>(0, n),
                     COMPOSE_LAMBDA (const Int& idx)
                     { q_min_1d(idx) = ko::max<Real>(q_min_1d(idx), 0); });
  }
  // Loop differently due to performance diff on CPU.
#ifdef COMPOSE_PORT
  const auto f = COMPOSE_LAMBDA (const Int& idx) {
    const Int ie = nets + idx/(nsuplev*qsize);
    const Int q = (idx / nsuplev) % qsize;
    const Int spli = idx % nsuplev;
#else
  for (Int ie = nets; ie <= nete; ++ie) {
#endif
    const auto spheremp1 = subview_ie(ie, spheremp);
    const auto dp3d_c1 = subview_ie(ie, dp3d_c);
    const auto qdp_p1 = subview_ie(ie, qdp_p);
    const auto q_c1 = subview_ie(ie, q_c);
#ifndef COMPOSE_PORT
    for (Int q = 0; q < qsize; ++q)
    for (Int spli = 0; spli < cdr.nsuplev; ++spli) {
#endif
    const Int k0 = nsublev*spli;
    const Int ti = cdr_over_super_levels ? q : spli*qsize + q;
    const bool nonneg = nonnegs[q];
    Real Qm = 0, Qm_min = 0, Qm_max = 0, Qm_prev = 0, rhom = 0, volume = 0;
    Int ie_idx;
    if (caas_in_suplev)
      ie_idx = (cdr_over_super_levels ?
                nsuplev*ie + spli :
                ie);
    for (Int sbli = 0; sbli < nsublev; ++sbli) {
      const auto k = k0 + sbli;
      if ( ! caas_in_suplev)
        ie_idx = (cdr_over_super_levels ?
                  nlevwrem*ie + k :
                  nsublev*ie + sbli);
      const auto lci = ie2lci[ie_idx];
      if ( ! caas_in_suplev) {
        Qm = 0; Qm_min = 0; Qm_max = 0; Qm_prev = 0;
        rhom = 0;
        volume = 0;
      }
      if (k < nlev) {
        for (Int g = 0; g < np2; ++g) {
          const auto smp = spheremp1(g);
          volume += smp;
          const Real rhomij = dp3d_c1(np1,g,k) * smp;
          rhom += rhomij;
          Qm += q_c1(q,g,k) * rhomij;
          auto& q_min_val = idx_qext(q_min,ie,q,g,k);
          if ( ! cedr::impl::OnGpu<typename MT::DES>::value && nonneg)
            q_min_val = ko::max<Real>(q_min_val, 0);
          Qm_min += q_min_val * rhomij;
          Qm_max += idx_qext(q_max,ie,q,g,k) * rhomij;
          Qm_prev += qdp_p1(n0_qdp,q,g,k) * smp;
        }
      }
      const bool write = ! caas_in_suplev || sbli == nsublev-1;
      if (write) {
        // For now, handle just one rhom. For feasible global problems,
        // it's used only as a weight vector in QLT, so it's fine. In fact,
        // use just the cell geometry, rather than total density, since in QLT
        // this field is used as a weight vector.
        //todo Generalize to one rhom field per level. Until then, we're not
        // getting QLT's safety benefit.
        if (ti == 0) cedr_cdr.set_rhom(lci, 0, volume);
        cedr_cdr.set_Qm(lci, ti, Qm, Qm_min, Qm_max, Qm_prev);
        if (Qm_prev < -0.5)
          warn_on_Qm_prev_negative<MT>(Qm_prev, rank, ie, ie2gci, np2, spli, k0, q,
                                       ti, sbli, lci, k, n0_qdp, np1, qdp_p, dp3d_c);
      }
    }
#ifdef COMPOSE_PORT
  };
  ko::fence();
  ko::parallel_for(ko::RangePolicy<typename MT::DES>(0, (nete - nets + 1)*nsuplev*qsize), f);
#else
  }}
#endif
}

template <typename MT>
void run_global (CDR<MT>& cdr, const Data& d, Real* q_min_r, const Real* q_max_r,
                 const Int nets, const Int nete) {
  if (dynamic_cast<typename CDR<MT>::QLTT*>(cdr.cdr.get()))
    run_global<4, MT, typename CDR<MT>::QLTT>(
      cdr, dynamic_cast<typename CDR<MT>::QLTT*>(cdr.cdr.get()),
      d, q_min_r, q_max_r, nets, nete);
  else if (dynamic_cast<typename CDR<MT>::CAAST*>(cdr.cdr.get()))
    run_global<4, MT, typename CDR<MT>::CAAST>(
      cdr, dynamic_cast<typename CDR<MT>::CAAST*>(cdr.cdr.get()),
      d, q_min_r, q_max_r, nets, nete);
  else
    cedr_throw_if(true, "run_global: could not cast cdr.");
  ko::fence();
  { Timer t("02_run_cdr");
    run_cdr(cdr); }
}

template void
run_global(CDR<ko::MachineTraits>& cdr, const Data& d, Real* q_min_r, const Real* q_max_r,
           const Int nets, const Int nete);

} // namespace sl
} // namespace homme
