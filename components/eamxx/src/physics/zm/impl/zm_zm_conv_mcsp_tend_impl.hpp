#ifndef ZM_ZM_CONV_MCSP_TEND_IMPL_HPP
#define ZM_ZM_CONV_MCSP_TEND_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_conv_mcsp_tend. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_conv_mcsp_tend(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Real& ztodt, // 2x physics time step
  const Int& jctop, // cloud top level indices
  const uview_1d<const Real>& state_pmid, // physics state mid-point pressure
  const uview_1d<const Real>& state_pint, // physics state interface pressure
  const uview_1d<const Real>& state_pdel, // physics state pressure thickness
  const uview_1d<const Real>& state_s, // physics state dry energy
  const uview_1d<const Real>& state_q, // physics state specific humidity
  const uview_1d<const Real>& state_u, // physics state u momentum
  const uview_1d<const Real>& state_v, // physics state v momentum
  const uview_1d<const Real>& ptend_zm_s, // input ZM tendency for dry energy (DSE)
  const uview_1d<const Real>& ptend_zm_q, // input ZM tendency for specific humidity (qv)
  // Inputs/Outputs
  const uview_1d<Real>& ptend_s, // output tendency of DSE
  const uview_1d<Real>& ptend_q, // output tendency of qv
  const uview_1d<Real>& ptend_u, // output tendency of u-wind
  const uview_1d<Real>& ptend_v, // output tendency of v-wind
  // Outputs
  const uview_1d<Real>& mcsp_dt_out, // final MCSP tendency for DSE
  const uview_1d<Real>& mcsp_dq_out, // final MCSP tendency for qv
  const uview_1d<Real>& mcsp_du_out, // final MCSP tendency for u wind
  const uview_1d<Real>& mcsp_dv_out, // final MCSP tendency for v wind
  Real& mcsp_freq, // MSCP frequency for output
  Real& mcsp_shear, // shear used to check against threshold
  Real& zm_depth) // pressure depth of ZM heating
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
