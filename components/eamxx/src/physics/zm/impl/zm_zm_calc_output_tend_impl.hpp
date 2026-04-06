#ifndef ZM_ZM_CALC_OUTPUT_TEND_IMPL_HPP
#define ZM_ZM_CALC_OUTPUT_TEND_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_calc_output_tend. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_calc_output_tend(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& jt, // level index of updraft top
  const Int& mx, // level index of updraft base
  const Real& dsubcld, // sub-cloud layer thickness
  const uview_1d<const Real>& p_del, // pressure thickness
  const uview_1d<const Real>& s_int, // ambient interface dry energy
  const uview_1d<const Real>& q_int, // ambient interface specific humidity
  const uview_1d<const Real>& s_upd, // updraft dry energy
  const uview_1d<const Real>& q_upd, // updraft specific humidity
  const uview_1d<const Real>& mflx_up, // updraft mass flux
  const uview_1d<const Real>& detr_up, // updraft detrainment
  const uview_1d<const Real>& mflx_dn, // downdraft mass flux
  const uview_1d<const Real>& s_dnd, // downdraft dry energy
  const uview_1d<const Real>& q_dnd, // downdraft specific humidity
  const uview_1d<const Real>& ql, // cloud liquid water
  const uview_1d<const Real>& evp, // evaporation
  const uview_1d<const Real>& cu, // updraft condensation
  // Outputs
  const uview_1d<Real>& dsdt, // output tendency for dry energy
  const uview_1d<Real>& dqdt, // output tendency for specific humidity
  const uview_1d<Real>& dl) // output tendency for cloud liquid water
{
  //----------------------------------------------------------------------------
  // find the highest level top and bottom levels of convection
  // trivial for single column: ktm = jt, kbm = mx
  const Int ktm = jt;
  const Int kbm = mx;

  //----------------------------------------------------------------------------
  // initialize variables
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg+1, pver), [&] (const Int& k) {
    dsdt(k) = 0;
    dqdt(k) = 0;
    dl(k)   = 0;
  });

  team.team_barrier();

  //----------------------------------------------------------------------------
  // calculate large-scale tendencies
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, ktm, pver-1), [&] (const Int& k) {
    const Real emc = -cu(k) + evp(k); // condensation in updraft and evaporating rain in downdraft

    dsdt(k) = -PC::LatVap.value/PC::Cpair.value*emc +
              ( +mflx_up(k+1)*(s_upd(k+1)-s_int(k+1)) - mflx_up(k)*(s_upd(k)-s_int(k))
                +mflx_dn(k+1)*(s_dnd(k+1)-s_int(k+1)) - mflx_dn(k)*(s_dnd(k)-s_int(k))
              )/p_del(k);

    dqdt(k) = emc +
              ( +mflx_up(k+1)*(q_upd(k+1)-q_int(k+1)) - mflx_up(k)*(q_upd(k)-q_int(k))
                +mflx_dn(k+1)*(q_dnd(k+1)-q_int(k+1)) - mflx_dn(k)*(q_dnd(k)-q_int(k))
              )/p_del(k);

    dl(k) = detr_up(k)*ql(k+1);
  });

  team.team_barrier();

  //----------------------------------------------------------------------------
  // calculate large-scale tendencies at and below cloud base
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    for (int k = kbm; k < pver; ++k) {
      if (k == mx) {
        dsdt(k) = (Real(1)/dsubcld)* ( -mflx_up(k)*(s_upd(k)-s_int(k))
                                       -mflx_dn(k)*(s_dnd(k)-s_int(k)) );
        dqdt(k) = (Real(1)/dsubcld)* ( -mflx_up(k)*(q_upd(k)-q_int(k))
                                       -mflx_dn(k)*(q_dnd(k)-q_int(k)) );
      } else if (k > mx) {
        dsdt(k) = dsdt(k-1);
        dqdt(k) = dqdt(k-1);
      }
    }
  });
}

} // namespace zm
} // namespace scream

#endif
