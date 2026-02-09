#ifndef ZM_ZM_TRANSPORT_TRACER_IMPL_HPP
#define ZM_ZM_TRANSPORT_TRACER_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_transport_tracer. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_transport_tracer(
  // Inputs
  const MemberType& team,
  const Int& pver,                        // number of mid-point levels
  const uview_1d<const bool>& doconvtran, // flag for doing convective transport
  const uview_2d<const Real>& q,          // tracer array (including water vapor)
  const Int& ncnst,                       // number of tracers to transport
  const uview_1d<const Real>& mu,         // mass flux up
  const uview_1d<const Real>& md,         // mass flux down
  const uview_1d<const Real>& du,         // mass detraining from updraft
  const uview_1d<const Real>& eu,         // mass entraining from updraft
  const uview_1d<const Real>& ed,         // mass entraining from downdraft
  const uview_1d<const Real>& dp,         // delta pressure between interfaces
  const Int& jt,                          // index of cloud top for each column
  const Int& mx,                          // index of cloud top for each column
  const Int& ideep,                       // gathering array
  const Int& il1g,                        // gathered min ncol index
  const Int& il2g,                        // gathered max ncol index
  const uview_2d<const Real>& fracis,     // fraction of tracer that is insoluble
  const uview_1d<const Real>& dpdry,      // delta pressure between interfaces
  const Real& dt,                         // model time increment)
  // Outputs
  const uview_2d<Real>& dqdt)             // output tendency array
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
