#include "compose_cedr_cdr.hpp"
#include "compose_cedr_sl.hpp"

namespace homme {
namespace sl {

template <typename MT>
void check (CDR<MT>& cdr, Data& d, const Real* q_min_r, const Real* q_max_r,
            const Int nets, const Int nete) {
  using cedr::mpi::reduce;

  const auto& ta = *d.ta;
  const Int np = ta.np, np2 = np*np, nlev = ta.nlev, nsuplev = cdr.nsuplev,
    qsize = ta.qsize, nprob = cdr.threed ? 1 : nsuplev;

  Kokkos::View<Real**, Kokkos::Serial>
    mass_p("mass_p", nprob, qsize), mass_c("mass_c", nprob, qsize),
    mass_lo("mass_lo", nprob, qsize), mass_hi("mass_hi", nprob, qsize),
    q_lo("q_lo", nprob, qsize), q_hi("q_hi", nprob, qsize),
    q_min_l("q_min_l", nprob, qsize), q_max_l("q_max_l", nprob, qsize),
    qd_lo("qd_lo", nprob, qsize), qd_hi("qd_hi", nprob, qsize);
  Kokkos::deep_copy(q_lo,  1e200);
  Kokkos::deep_copy(q_hi, -1e200);
  Kokkos::deep_copy(q_min_l,  1e200);
  Kokkos::deep_copy(q_max_l, -1e200);
  Kokkos::deep_copy(qd_lo, 0);
  Kokkos::deep_copy(qd_hi, 0);

#ifdef COMPOSE_PORT
  const auto q_min = ko::create_mirror_view(ta.q_min);
  const auto q_max = ko::create_mirror_view(ta.q_max);
  ko::deep_copy(q_min, ta.q_min);
  ko::deep_copy(q_max, ta.q_max);
#else
  const QExtremaHConst<ko::MachineTraits>
    q_min(q_min_r, ta.nelemd, ta.qsize, ta.nlev, ta.np2),
    q_max(q_max_r, ta.nelemd, ta.qsize, ta.nlev, ta.np2);
#endif
  const auto np1 = ta.np1;
  const auto n0_qdp = ta.n0_qdp;
  const auto n1_qdp = ta.n1_qdp;
#ifdef COMPOSE_PORT
  const auto spheremp = ko::cmvdc(ta.spheremp);
  const auto dp3d_c = ko::cmvdc(ta.dp3d);
  const auto qdp_pc = ko::cmvdc(ta.qdp);
  const auto q_c = ko::cmvdc(ta.q);
#else
  const auto& spheremp = ta.pspheremp;
  const auto& dp3d_c = ta.pdp3d;
  const auto& qdp_pc = ta.pqdp;
  const auto& q_c = ta.pq;
#endif

  Int iprob = 0;

  bool fp_issue = false; // Limit output once the first issue is seen.
  for (Int ie = nets; ie <= nete; ++ie) {
    for (Int spli = 0; spli < nsuplev; ++spli) {
      if (nprob > 1) iprob = spli;
      for (Int k = spli*cdr.nsublev; k < (spli+1)*cdr.nsublev; ++k) {
        if (k >= nlev) continue;
        if ( ! fp_issue) {
          for (Int g = 0; g < np2; ++g) {
            // FP issues.
            if (std::isnan(dp3d_c(ie,np1,g,k)))
            { pr("dp3d NaN:" pu(k) pu(g)); fp_issue = true; }
            if (std::isinf(dp3d_c(ie,np1,g,k)))
            { pr("dp3d Inf:" pu(k) pu(g)); fp_issue = true; }
          }
        }
        for (Int q = 0; q < qsize; ++q) {
          Real qlo_s = idx_qext(q_min,ie,q,0,k), qhi_s = idx_qext(q_max,ie,q,0,k);
          for (Int g = 0; g < np2; ++g) {
            qlo_s = std::min(qlo_s, idx_qext(q_min,ie,q,g,k));
            qhi_s = std::max(qhi_s, idx_qext(q_max,ie,q,g,k));
          }
          for (Int g = 0; g < np2; ++g) {
            // FP issues.
            if ( ! fp_issue) {
              for (Int i_qdp : {0, 1}) {
                const Int n_qdp = i_qdp == 0 ? n0_qdp : n1_qdp;
                if (std::isnan(qdp_pc(ie,n_qdp,q,g,k)))
                { pr("qdp NaN:" puf(i_qdp) pu(q) pu(k) pu(g)); fp_issue = true; }
                if (std::isinf(qdp_pc(ie,n_qdp,q,g,k)))
                { pr("qdp Inf:" puf(i_qdp) pu(q) pu(k) pu(g)); fp_issue = true; }
              }
              if (std::isnan(q_c(ie,q,g,k)))
              { pr("q NaN:" pu(q) pu(k) pu(g)); fp_issue = true; }
              if (std::isinf(q_c(ie,q,g,k)))
              { pr("q Inf:" pu(q) pu(k) pu(g)); fp_issue = true; }
            }
            // Mass conservation.
            mass_p(iprob,q) += qdp_pc(ie,n0_qdp,q,g,k) * spheremp(ie,g);
            mass_c(iprob,q) += qdp_pc(ie,n1_qdp,q,g,k) * spheremp(ie,g);
            // Local bound constraints w.r.t. cell-local extrema.
            if (q_c(ie,q,g,k) < qlo_s)
              qd_lo(iprob,q) = std::max(qd_lo(iprob,q), qlo_s - q_c(ie,q,g,k));
            if (q_c(ie,q,g,k) > qhi_s)
              qd_hi(iprob,q) = std::max(qd_hi(iprob,q), q_c(ie,q,g,k) - qhi_s);
            // Safety problem bound constraints.
            mass_lo(iprob,q) += (idx_qext(q_min,ie,q,g,k) * dp3d_c(ie,np1,g,k) *
                                 spheremp(ie,g));
            mass_hi(iprob,q) += (idx_qext(q_max,ie,q,g,k) * dp3d_c(ie,np1,g,k) *
                                 spheremp(ie,g));
            q_lo(iprob,q) = std::min(q_lo(iprob,q), q_c(ie,q,g,k));
            q_hi(iprob,q) = std::max(q_hi(iprob,q), q_c(ie,q,g,k));
            q_min_l(iprob,q) = std::min(q_min_l(iprob,q), idx_qext(q_min,ie,q,g,k));
            q_max_l(iprob,q) = std::max(q_max_l(iprob,q), idx_qext(q_max,ie,q,g,k));
          }
        }
      }
    }
  }

#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    if ( ! d.check)
      d.check = std::make_shared<Data::Check>(nprob, qsize);
    auto& c = *d.check;
    Kokkos::deep_copy(c.mass_p, 0);
    Kokkos::deep_copy(c.mass_c, 0);
    Kokkos::deep_copy(c.mass_lo, 0);
    Kokkos::deep_copy(c.mass_hi, 0);
    Kokkos::deep_copy(c.q_lo,  1e200);
    Kokkos::deep_copy(c.q_hi, -1e200);
    Kokkos::deep_copy(c.q_min_l,  1e200);
    Kokkos::deep_copy(c.q_max_l, -1e200);
    Kokkos::deep_copy(c.qd_lo, 0);
    Kokkos::deep_copy(c.qd_hi, 0);
  }

#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp critical
#endif
  {
    auto& c = *d.check;
    for (Int spli = 0; spli < nprob; ++spli) {
      if (nprob > 1) iprob = spli;
      for (Int q = 0; q < qsize; ++q) {
        c.mass_p(iprob,q) += mass_p(iprob,q);
        c.mass_c(iprob,q) += mass_c(iprob,q);
        c.qd_lo(iprob,q) = std::max(c.qd_lo(iprob,q), qd_lo(iprob,q));
        c.qd_hi(iprob,q) = std::max(c.qd_hi(iprob,q), qd_hi(iprob,q));
        c.mass_lo(iprob,q) += mass_lo(iprob,q);
        c.mass_hi(iprob,q) += mass_hi(iprob,q);
        c.q_lo(iprob,q) = std::min(c.q_lo(iprob,q), q_lo(iprob,q));
        c.q_hi(iprob,q) = std::max(c.q_hi(iprob,q), q_hi(iprob,q));
        c.q_min_l(iprob,q) = std::min(c.q_min_l(iprob,q), q_min_l(iprob,q));
        c.q_max_l(iprob,q) = std::max(c.q_max_l(iprob,q), q_max_l(iprob,q));
      }
    }
  }

#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp barrier
# pragma omp master
#endif
  {
    Kokkos::View<Real**, Kokkos::Serial>
      mass_p_g("mass_p_g", nprob, qsize), mass_c_g("mass_c_g", nprob, qsize),
      mass_lo_g("mass_lo_g", nprob, qsize), mass_hi_g("mass_hi_g", nprob, qsize),
      q_lo_g("q_lo_g", nprob, qsize), q_hi_g("q_hi_g", nprob, qsize),
      q_min_g("q_min_g", nprob, qsize), q_max_g("q_max_g", nprob, qsize),
      qd_lo_g("qd_lo_g", nprob, qsize), qd_hi_g("qd_hi_g", nprob, qsize);

    const auto& p = *cdr.p;
    const auto& c = *d.check;
    const auto root = cdr.p->root();
    const auto N = nprob*qsize;

    reduce(p, c.mass_p.data(), mass_p_g.data(), N, MPI_SUM, root);
    reduce(p, c.mass_c.data(), mass_c_g.data(), N, MPI_SUM, root);
    reduce(p, c.qd_lo.data(), qd_lo_g.data(), N, MPI_MAX, root);
    reduce(p, c.qd_hi.data(), qd_hi_g.data(), N, MPI_MAX, root);
    // Safety problem.
    reduce(p, c.mass_lo.data(), mass_lo_g.data(), N, MPI_SUM, root);
    reduce(p, c.mass_hi.data(), mass_hi_g.data(), N, MPI_SUM, root);
    reduce(p, c.q_lo.data(), q_lo_g.data(), N, MPI_MIN, root);
    reduce(p, c.q_hi.data(), q_hi_g.data(), N, MPI_MAX, root);
    reduce(p, c.q_min_l.data(), q_min_g.data(), N, MPI_MIN, root);
    reduce(p, c.q_max_l.data(), q_max_g.data(), N, MPI_MAX, root);

    if (cdr.p->amroot()) {
      const Real tol = 1e4*std::numeric_limits<Real>::epsilon();
      for (Int k = 0; k < nprob; ++k)
        for (Int q = 0; q < qsize; ++q) {
          const Real rd = cedr::util::reldif(mass_p_g(k,q), mass_c_g(k,q));
          if (rd > tol)
            pr(puf(k) pu(q) pu(mass_p_g(k,q)) pu(mass_c_g(k,q)) pu(rd));
          if (mass_lo_g(k,q) <= mass_c_g(k,q) && mass_c_g(k,q) <= mass_hi_g(k,q)) {
            // Local problems should be feasible.
            if (qd_lo_g(k,q) > 0)
              pr(puf(k) pu(q) pu(qd_lo_g(k,q)));
            if (qd_hi_g(k,q) > 0)
              pr(puf(k) pu(q) pu(qd_hi_g(k,q)));
          } else {
            // Safety problem must hold.
            if (q_lo_g(k,q) < q_min_g(k,q))
              pr(puf(k) pu(q) pu(q_lo_g(k,q) - q_min_g(k,q)) pu(q_min_g(k,q)));
            if (q_hi_g(k,q) > q_max_g(k,q))
              pr(puf(k) pu(q) pu(q_max_g(k,q) - q_hi_g(k,q)) pu(q_max_g(k,q)));
          }
        }
    }
  }
}

template void
check(CDR<ko::MachineTraits>& cdr, Data& d, const Real* q_min_r, const Real* q_max_r,
      const Int nets, const Int nete);

} // namespace sl
} // namespace homme
