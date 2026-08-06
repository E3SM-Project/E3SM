#include "compose_cedr_qlt.hpp"
#include "compose_cedr.hpp"
#include "compose_kokkos.hpp"

#define THREAD_QLT_RUN

namespace ko = Kokkos;

namespace homme {
namespace compose {

typedef cedr::Int Int;
typedef cedr::Real Real;
template <typename ES> using RealList = typename cedr::qlt::QLT<ES>::RealList;

//todo All of this VerticalLevelsData-related code should be impl'ed
// using a new optional root-node-function QLT registration function.
struct VerticalLevelsData {
  typedef std::shared_ptr<VerticalLevelsData> Ptr;

  RealList<ko::Serial> lo, hi, mass, ones, wrk;

  VerticalLevelsData (const cedr::Int n)
    : lo("lo", n), hi("hi", n), mass("mass", n), ones("ones", n), wrk("wrk", n)
  {
    for (cedr::Int k = 0; k < n; ++k) ones(k) = 1;
  }
};

static Int solve (const Int n, const Real* a, const Real b,
                  const Real* xlo, const Real* xhi,
                  Real* x, Real* wrk) {
#ifndef NDEBUG
  cedr_assert(b >= 0);
  for (Int i = 0; i < n; ++i) {
    cedr_assert(a[i] > 0);
    cedr_assert(xlo[i] >= 0);
    cedr_assert(xhi[i] >= xlo[i]);
  }
#endif
  Int status = 0;
  Real tot_lo = 0, tot_hi = 0;
  for (Int i = 0; i < n; ++i) tot_lo += a[i]*xlo[i];
  for (Int i = 0; i < n; ++i) tot_hi += a[i]*xhi[i];
  if (b < tot_lo) {
    status = -2;
    for (Int i = 0; i < n; ++i) wrk[i] = 0;
    for (Int i = 0; i < n; ++i) x[i] = xlo[i];
    // Find a new xlo >= 0 minimally far from the current one. This
    // is also the solution x.
    cedr::local::caas(n, a, b, wrk, xlo, x, x, false);
  } else if (b > tot_hi) {
    status = -1;
    const Real f = b/tot_hi;
    // a[i] divides out.
    for (Int i = 0; i < n; ++i) x[i] = f*xhi[i];
  } else {
    cedr::local::caas(n, a, b, xlo, xhi, x, x, false);
  }
  return status;
}

static Int solve (const Int nlev, const VerticalLevelsData& vld,
                  const Real& tot_mass) {
  return solve(nlev, vld.ones.data(), tot_mass, vld.lo.data(), vld.hi.data(),
               vld.mass.data(), vld.wrk.data());    
}

static Int solve_unittest () {
  static const auto eps = std::numeric_limits<Real>::epsilon();
  static const Int n = 7;

  Real a[n], xlo[n], xhi[n], x[n], wrk[n];
  static const Real x0  [n] = { 1.2, 0.5,3  , 2  , 1.5, 1.8,0.2};
  static const Real dxlo[n] = {-0.1,-0.2,0.5,-1.5,-0.1,-1.1,0.1};
  static const Real dxhi[n] = { 0.1,-0.1,1  ,-0.5, 0.1,-0.2,0.5};
  for (Int i = 0; i < n; ++i) a[i] = i+1;
  for (Int i = 0; i < n; ++i) xlo[i] = x0[i] + dxlo[i];
  for (Int i = 0; i < n; ++i) xhi[i] = x0[i] + dxhi[i];
  Real b, b1;
  Int status, nerr = 0;

  const auto check_mass = [&] () {
    b1 = 0;
    for (Int i = 0; i < n; ++i) b1 += a[i]*x[i];
    if (std::abs(b1 - b) >= 10*eps*b) ++nerr;
  };

  for (Int i = 0; i < n; ++i) x[i] = x0[i];
  b = 0;
  for (Int i = 0; i < n; ++i) b += a[i]*xlo[i];
  b *= 0.9;
  status = solve(n, a, b, xlo, xhi, x, wrk);
  if (status != -2) ++nerr;
  check_mass();
  for (Int i = 0; i < n; ++i) if (x[i] > xhi[i]*(1 + 10*eps)) ++nerr;

  for (Int i = 0; i < n; ++i) x[i] = x0[i];
  b = 0;
  for (Int i = 0; i < n; ++i) b += a[i]*xhi[i];
  b *= 1.1;
  status = solve(n, a, b, xlo, xhi, x, wrk);
  if (status != -1) ++nerr;
  check_mass();
  for (Int i = 0; i < n; ++i) if (x[i] < xlo[i]*(1 - 10*eps)) ++nerr;

  for (Int i = 0; i < n; ++i) x[i] = x0[i];
  b = 0;
  for (Int i = 0; i < n; ++i) b += 0.5*a[i]*(xlo[i] + xhi[i]);
  status = solve(n, a, b, xlo, xhi, x, wrk);
  if (status != 0) ++nerr;
  check_mass();
  for (Int i = 0; i < n; ++i) if (x[i] < xlo[i]*(1 - 10*eps)) ++nerr;
  for (Int i = 0; i < n; ++i) if (x[i] > xhi[i]*(1 + 10*eps)) ++nerr;

  return nerr;
}

template <typename ES>
void QLT<ES>::reconcile_vertical (const Int problem_type, const Int bd_os,
                                  const Int bis, const Int bie) {
  using cedr::ProblemType;

  cedr_assert((problem_type & ProblemType::shapepreserve) &&
              (problem_type & ProblemType::conserve));

  auto& md = this->o.md_;
  auto& bd = this->o.bd_;
  const auto& vld = *vld_;
  const Int nlev = vld.lo.extent_int(0);
  const Int nprob = (bie - bis)/nlev;

#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
# pragma omp master
#endif
  for (Int pi = 0; pi < nprob; ++pi) {
    const Int bd_os_pi = bd_os + md.a_d.trcr2bl2r(md.a_d.bidx2trcr(bis + pi));
#ifdef RV_DIAG
    Real oob = 0, tot_mass_slv = 0, oob_slv = 0;
#endif
    Real tot_mass = 0;
    for (Int k = 0; k < nlev; ++k) {
      const Int bd_os_k = bd_os_pi + nprob*4*k;
      vld.lo  (k) = bd.l2r_data(bd_os_k    );
      vld.hi  (k) = bd.l2r_data(bd_os_k + 2);
      vld.mass(k) = bd.l2r_data(bd_os_k + 3); // previous mass, not current one
      tot_mass += vld.mass(k);
#ifdef RV_DIAG
      if (vld.mass(k) < vld.lo(k)) oob += vld.lo(k) - vld.mass(k);
      if (vld.mass(k) > vld.hi(k)) oob += vld.mass(k) - vld.hi(k);
#endif
    }
    solve(nlev, vld, tot_mass);
    for (Int k = 0; k < nlev; ++k) {
      const Int bd_os_k = bd_os_pi + nprob*4*k;
      bd.l2r_data(bd_os_k + 3) = vld.mass(k); // previous mass, not current one
#ifdef RV_DIAG
      tot_mass_slv += vld.mass(k);
      if (vld.mass(k) < vld.lo(k)) oob_slv += vld.lo(k) - vld.mass(k);
      if (vld.mass(k) > vld.hi(k)) oob_slv += vld.mass(k) - vld.hi(k);
#endif
    }
#ifdef RV_DIAG
    printf("%2d %9.2e %9.2e %9.2e %9.2e\n", pi,
           oob/tot_mass, oob_slv/tot_mass,
           tot_mass, std::abs(tot_mass_slv - tot_mass)/tot_mass);
#endif
  }
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
# pragma omp barrier
#endif
}

template <typename ES>
QLT<ES>::QLT (const cedr::mpi::Parallel::Ptr& p, const cedr::Int& ncells,
              const cedr::tree::Node::Ptr& tree, const cedr::CDR::Options& options,
              const cedr::Int& vertical_levels)
  : cedr::qlt::QLT<ES>(p, ncells, tree, options)
{
  if (vertical_levels) {
    if (ko::OnGpu<ES>::value)
      cedr_throw_if(ko::OnGpu<ES>::value,
                    "QLT does not yet support vertical_levels on gpu.");
    else
      vld_ = std::make_shared<VerticalLevelsData>(vertical_levels);
  }
}

template <typename ES>
void QLT<ES>::run () {
  if (ko::OnGpu<ES>::value)
    Super::run();
  else
    runimpl();
}

template <typename ES>
void QLT<ES>::runimpl () {
  static const int mpitag = 42;
  using cedr::Int;
  using cedr::Real;
  using cedr::ProblemType;
  using cedr::tree::NodeSets;
  namespace mpi = cedr::mpi;
  auto& md_ = this->o.md_;
  auto& bd_ = this->o.bd_;
  auto& ns_ = this->ns_;
  auto& p_ = this->p_;
#if ! defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
# pragma omp master
  {
#endif
    // Number of data per slot.
    const Int l2rndps = md_.a_d.prob2bl2r[md_.nprobtypes];
    const Int r2lndps = md_.a_d.prob2br2l[md_.nprobtypes];

    // Leaves to root.
    for (size_t il = 0; il < ns_->levels.size(); ++il) {
      auto& lvl = ns_->levels[il];

      // Set up receives.
      if (lvl.kids.size()) {
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#       pragma omp master
#endif
        {
          for (size_t i = 0; i < lvl.kids.size(); ++i) {
            const auto& mmd = lvl.kids[i];
            mpi::irecv(*p_, &bd_.l2r_data(mmd.offset*l2rndps), mmd.size*l2rndps, mmd.rank,
                       mpitag, &lvl.kids_req[i]);
          }
          mpi::waitall(lvl.kids_req.size(), lvl.kids_req.data());
        }
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#       pragma omp barrier
#endif
      }

      // Combine kids' data.
      if (lvl.nodes.size()) {
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#       pragma omp for
#endif
        for (size_t ni = 0; ni < lvl.nodes.size(); ++ni) {
          const auto lvlidx = lvl.nodes[ni];
          const auto n = ns_->node_h(lvlidx);
          if ( ! n->nkids) continue;
          cedr_kernel_assert(n->nkids == 2);
          // Total density.
          bd_.l2r_data(n->offset*l2rndps) =
            (bd_.l2r_data(ns_->node_h(n->kids[0])->offset*l2rndps) +
             bd_.l2r_data(ns_->node_h(n->kids[1])->offset*l2rndps));
          // Tracers.
          for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
            const Int problem_type = md_.get_problem_type(pti);
            const bool nonnegative = problem_type & ProblemType::nonnegative;
            const bool shapepreserve = problem_type & ProblemType::shapepreserve;
            const bool conserve = problem_type & ProblemType::conserve;
            const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
#if defined THREAD_QLT_RUN && defined COMPOSE_COLUMN_OPENMP
#           pragma omp parallel for
#endif
            for (Int bi = bis; bi < bie; ++bi) {
              const Int bdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
              Real* const me = &bd_.l2r_data(n->offset*l2rndps + bdi);
              const auto kid0 = ns_->node_h(n->kids[0]);
              const auto kid1 = ns_->node_h(n->kids[1]);
              const Real* const k0 = &bd_.l2r_data(kid0->offset*l2rndps + bdi);
              const Real* const k1 = &bd_.l2r_data(kid1->offset*l2rndps + bdi);
              if (nonnegative) {
                me[0] = k0[0] + k1[0];
                if (conserve) me[1] = k0[1] + k1[1];
              } else {
                me[0] = shapepreserve ? k0[0] + k1[0] : cedr::impl::min(k0[0], k1[0]);
                me[1] = k0[1] + k1[1];
                me[2] = shapepreserve ? k0[2] + k1[2] : cedr::impl::max(k0[2], k1[2]);
                if (conserve) me[3] = k0[3] + k1[3] ;
              }
            }
          }
        }
      }

      // Send to parents.
      if (lvl.me.size())
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#     pragma omp master
#endif
      {
        for (size_t i = 0; i < lvl.me.size(); ++i) {
          const auto& mmd = lvl.me[i];
          mpi::isend(*p_, &bd_.l2r_data(mmd.offset*l2rndps), mmd.size*l2rndps,
                     mmd.rank, mpitag);
        }
      }

#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#     pragma omp barrier
#endif
    }

    // Root.
    if ( ! (ns_->levels.empty() || ns_->levels.back().nodes.size() != 1 ||
            ns_->node_h(ns_->levels.back().nodes[0])->parent >= 0)) {
      const auto n = ns_->node_h(ns_->levels.back().nodes[0]);
      for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
        const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
        if (bie == bis) continue;
        const Int problem_type = md_.get_problem_type(pti);
        if (vld_) reconcile_vertical(problem_type, n->offset*l2rndps, bis, bie);
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP && defined COMPOSE_COLUMN_OPENMP
#       pragma omp parallel
#endif
#if defined THREAD_QLT_RUN && (defined COMPOSE_HORIZ_OPENMP || defined COMPOSE_COLUMN_OPENMP)
#       pragma omp for
#endif
        for (Int bi = bis; bi < bie; ++bi) {
          const Int l2rbdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
          const Int r2lbdi = md_.a_d.trcr2br2l(md_.a_d.bidx2trcr(bi));
          // If QLT is enforcing global mass conservation, set the root's r2l Qm
          // value to the l2r Qm_prev's sum; otherwise, copy the l2r Qm value to
          // the r2l one.
          const Int os = (problem_type & ProblemType::conserve ?
                          Super::MetaData::get_problem_type_l2r_bulk_size(problem_type) - 1 :
                          (problem_type & ProblemType::nonnegative ? 0 : 1));
          bd_.r2l_data(n->offset*r2lndps + r2lbdi) =
            bd_.l2r_data(n->offset*l2rndps + l2rbdi + os);
          if ((problem_type & ProblemType::consistent) &&
              ! (problem_type & ProblemType::shapepreserve)) {
            // Consistent but not shape preserving, so we're solving a dynamic range
            // preservation problem. We now know the global q_{min,max}. Start
            // propagating it leafward.
            bd_.r2l_data(n->offset*r2lndps + r2lbdi + 1) =
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 0);
            bd_.r2l_data(n->offset*r2lndps + r2lbdi + 2) =
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 2);
          }
        }
      }
    }

    // Root to leaves.
    for (size_t il = ns_->levels.size(); il > 0; --il) {
      auto& lvl = ns_->levels[il-1];

      if (lvl.me.size()) {
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#       pragma omp master
#endif
        {
          for (size_t i = 0; i < lvl.me.size(); ++i) {
            const auto& mmd = lvl.me[i];
            mpi::irecv(*p_, &bd_.r2l_data(mmd.offset*r2lndps), mmd.size*r2lndps, mmd.rank,
                       mpitag, &lvl.me_recv_req[i]);
          }
          mpi::waitall(lvl.me_recv_req.size(), lvl.me_recv_req.data());
        }
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#       pragma omp barrier
#endif
      }

      // Solve QP for kids' values.
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#     pragma omp for
#endif
      for (size_t ni = 0; ni < lvl.nodes.size(); ++ni) {
        const auto lvlidx = lvl.nodes[ni];
        const auto n = ns_->node_h(lvlidx);
        if ( ! n->nkids) continue;
        for (Int pti = 0; pti < md_.nprobtypes; ++pti) {
          const Int problem_type = md_.get_problem_type(pti);
          const Int bis = md_.a_d.prob2trcrptr[pti], bie = md_.a_d.prob2trcrptr[pti+1];
#if defined THREAD_QLT_RUN && defined COMPOSE_COLUMN_OPENMP
#         pragma omp parallel for
#endif
          for (Int bi = bis; bi < bie; ++bi) {
            const Int l2rbdi = md_.a_d.trcr2bl2r(md_.a_d.bidx2trcr(bi));
            const Int r2lbdi = md_.a_d.trcr2br2l(md_.a_d.bidx2trcr(bi));
            cedr_assert(n->nkids == 2);
            if ((problem_type & ProblemType::consistent) &&
                ! (problem_type & ProblemType::shapepreserve)) {
              // Pass q_{min,max} info along. l2r data are updated for use in
              // solve_node_problem. r2l data are updated for use in isend.
              const Real q_min = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 1);
              const Real q_max = bd_.r2l_data(n->offset*r2lndps + r2lbdi + 2);
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 0) = q_min;
              bd_.l2r_data(n->offset*l2rndps + l2rbdi + 2) = q_max;
              for (Int k = 0; k < 2; ++k) {
                const auto os = ns_->node_h(n->kids[k])->offset;
                bd_.l2r_data(os*l2rndps + l2rbdi + 0) = q_min;
                bd_.l2r_data(os*l2rndps + l2rbdi + 2) = q_max;
                bd_.r2l_data(os*r2lndps + r2lbdi + 1) = q_min;
                bd_.r2l_data(os*r2lndps + r2lbdi + 2) = q_max;
              }
            }
            const auto k0 = ns_->node_h(n->kids[0]);
            const auto k1 = ns_->node_h(n->kids[1]);
            cedr::qlt::impl::solve_node_problem(
              problem_type,
               bd_.l2r_data( n->offset*l2rndps),
              &bd_.l2r_data( n->offset*l2rndps + l2rbdi),
               bd_.r2l_data( n->offset*r2lndps + r2lbdi),
               bd_.l2r_data(k0->offset*l2rndps),
              &bd_.l2r_data(k0->offset*l2rndps + l2rbdi),
               bd_.r2l_data(k0->offset*r2lndps + r2lbdi),
               bd_.l2r_data(k1->offset*l2rndps),
              &bd_.l2r_data(k1->offset*l2rndps + l2rbdi),
               bd_.r2l_data(k1->offset*r2lndps + r2lbdi),
              this->options_.prefer_numerical_mass_conservation_to_numerical_bounds);
          }
        }
      }

      // Send.
      if (lvl.kids.size())
#if defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
#     pragma omp master
#endif
      {
        for (size_t i = 0; i < lvl.kids.size(); ++i) {
          const auto& mmd = lvl.kids[i];
          mpi::isend(*p_, &bd_.r2l_data(mmd.offset*r2lndps), mmd.size*r2lndps,
                     mmd.rank, mpitag);
        }
      }
    }
#if ! defined THREAD_QLT_RUN && defined COMPOSE_HORIZ_OPENMP
  }
#endif
}

template <typename ES>
Int QLT<ES>::unittest () {
  return solve_unittest();
}

template class QLT<ko::MachineTraits::DES>;

} // namespace compose
} // namespace homme
