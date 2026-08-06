// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_CAAS_INL_HPP
#define INCLUDE_CEDR_CAAS_INL_HPP

#include "cedr_util.hpp"

namespace cedr {
// ClipAndAssuredSum.
namespace caas {

template <typename ES> KOKKOS_INLINE_FUNCTION
void CAAS<ES>::DeviceOp
::set_rhom (const Int& lclcellidx, const Int& rhomidx, const Real& rhom) const {
  cedr_kernel_assert(lclcellidx >= 0 && lclcellidx < nlclcells_);
  cedr_kernel_assert(rhomidx >= 0 && rhomidx < nrhomidxs_);
  d_(lclcellidx) = rhom;
}

template <typename ES> KOKKOS_INLINE_FUNCTION
void CAAS<ES>::DeviceOp
::set_Qm (const Int& lclcellidx, const Int& tracer_idx,
          const Real& Qm, const Real& Qm_min, const Real& Qm_max,
          const Real Qm_prev) const {
  cedr_kernel_assert(lclcellidx >= 0 && lclcellidx < nlclcells_);
  cedr_kernel_assert(tracer_idx >= 0 && tracer_idx < probs_.extent_int(0));
  const Int nt = probs_.size();
  d_((1 +        tracer_idx)*nlclcells_ + lclcellidx) = Qm;
  d_((1 +   nt + tracer_idx)*nlclcells_ + lclcellidx) = Qm_min;
  d_((1 + 2*nt + tracer_idx)*nlclcells_ + lclcellidx) = Qm_max;
  if (need_conserve_)
    d_((1 + 3*nt + tracer_idx)*nlclcells_ + lclcellidx) = Qm_prev;
}

template <typename ES> KOKKOS_INLINE_FUNCTION
Real CAAS<ES>::DeviceOp::
get_Qm (const Int& lclcellidx, const Int& tracer_idx) const {
  cedr_kernel_assert(lclcellidx >= 0 && lclcellidx < nlclcells_);
  cedr_kernel_assert(tracer_idx >= 0 && tracer_idx < probs_.extent_int(0));
  return d_((1 + tracer_idx)*nlclcells_ + lclcellidx);
}

template <typename RealList, typename IntList>
KOKKOS_INLINE_FUNCTION static void
calc_Qm_scalars (const RealList& d, const IntList& probs,
                 const Int& nt, const Int& nlclcells,
                 const Int& k, const Int& os, const Int& i,
                 Real& Qm_clip, Real& Qm_term) {
  const Real Qm = d(os+i);
  Qm_term = (probs(k) & ProblemType::conserve ?
             d(os + nlclcells*3*nt + i) /* Qm_prev */ :
             Qm);
  const Real Qm_min = d(os + nlclcells*  nt + i);
  const Real Qm_max = d(os + nlclcells*2*nt + i);
  Qm_clip = cedr::impl::min(Qm_max, cedr::impl::max(Qm_min, Qm));
}

} // namespace caas
} // namespace cedr

#endif
