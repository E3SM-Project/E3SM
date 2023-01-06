#ifndef COMPOSE_PORT
#include "compose_cedr_caas.hpp"
#include "compose_kokkos.hpp"

namespace homme {
namespace compose {

typedef cedr::Int Int;
typedef cedr::Real Real;

// On Chrysalis, optimization level above O1 leads to BFB variance w.r.t. PE
// layout on the WC v2 problem, despite the e3sm_atm_developer and
// e3sm_atm_integration passing against baselines.
#ifdef __INTEL_COMPILER
# pragma intel optimization_level 1
#endif
void CAAS::run_horiz_omp () {
  cedr_assert(finished_setup_);
  cedr_assert(user_reducer_ != nullptr);
  reduce_locally_horiz_omp();
  (*user_reducer_)(*p_, send_.data(), recv_.data(),
                   o.nlclcells_ / user_reducer_->n_accum_in_place(),
                   recv_.size(), MPI_SUM);
  finish_locally_horiz_omp();
}

template <typename Func>
void homme_parallel_for (const Int& beg, const Int& end, const Func& f) {
#ifdef COMPOSE_HORIZ_OPENMP
# pragma omp for
#endif
  for (Int i = beg; i < end; ++i)
    f(i);
}

void CAAS::reduce_locally_horiz_omp () {
  const bool user_reduces = user_reducer_ != nullptr;
  cedr_assert(user_reduces); // assumption in Homme
  ConstExceptGnu Int nt = o.probs_.size(), nlclcells = o.nlclcells_;

  const auto& probs = o.probs_;
  const auto& send = send_;
  const auto& d = o.d_;
  const Int n_accum_in_place = user_reducer_->n_accum_in_place();
  const Int nlclaccum = nlclcells / n_accum_in_place;
  const auto calc_Qm_clip = COMPOSE_LAMBDA (const Int& j) {
    const auto k = j / nlclaccum;
    const auto bi = j % nlclaccum;
    const auto os = (k+1)*nlclcells;
    Real accum_clip = 0, accum_term = 0;
    for (Int ai = 0; ai < n_accum_in_place; ++ai) {
      const Int i = n_accum_in_place*bi + ai;
      Real Qm_clip, Qm_term;
      cedr::caas::calc_Qm_scalars(d, probs, nt, nlclcells, k, os, i, Qm_clip,
                                  Qm_term);
      d(os + i) = Qm_clip;
      accum_clip += Qm_clip;
      accum_term += Qm_term;
    }
    send(nlclaccum*      k  + bi) = accum_clip;
    send(nlclaccum*(nt + k) + bi) = accum_term;
  };
  homme_parallel_for(0, nt*nlclaccum, calc_Qm_clip);
  const auto set_Qm_minmax = COMPOSE_LAMBDA (const Int& j) {
    const auto k = 2*nt + j / nlclaccum;
    const auto bi = j % nlclaccum;
    const auto os = (k-nt+1)*nlclcells;
    Real accum_ext = 0;
    for (Int ai = 0; ai < n_accum_in_place; ++ai) {
      const Int i = n_accum_in_place*bi + ai;
      accum_ext += d(os + i);
    }
    send(nlclaccum*k + bi) = accum_ext;
  };
  homme_parallel_for(0, 2*nt*nlclaccum, set_Qm_minmax);
}

void CAAS::finish_locally_horiz_omp () {
  ConstExceptGnu Int nt = o.probs_.size(), nlclcells = o.nlclcells_;
  const auto& recv = recv_;
  const auto& d = o.d_;
  const auto adjust_Qm = COMPOSE_LAMBDA (const Int& k) {
    const auto os = (k+1)*nlclcells;
    const auto Qm_clip_sum = recv(     k);
    const auto Qm_sum      = recv(nt + k);
    const auto m = Qm_sum - Qm_clip_sum;
    if (m < 0) {
      const auto Qm_min_sum = recv(2*nt + k);
      auto fac = Qm_clip_sum - Qm_min_sum;
      if (fac > 0) {
        fac = m/fac;
        for (Int i = 0; i < nlclcells; ++i) {
          const auto Qm_min = d(os + nlclcells * nt + i);
          auto& Qm = d(os+i);
          Qm += fac*(Qm - Qm_min);
          Qm = cedr::impl::max(Qm_min, Qm);
        };
      }
    } else if (m > 0) {
      const auto Qm_max_sum = recv(3*nt + k);
      auto fac = Qm_max_sum - Qm_clip_sum;
      if (fac > 0) {
        fac = m/fac;
        for (Int i = 0; i < nlclcells; ++i) {
          const auto Qm_max = d(os + nlclcells*2*nt + i);
          auto& Qm = d(os+i);
          Qm += fac*(Qm_max - Qm);
          Qm = cedr::impl::min(Qm_max, Qm);
        };
      }
    }
  };
  homme_parallel_for(0, nt, adjust_Qm);  
}

} // namespace compose
} // namespace homme
#endif
