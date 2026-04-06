#ifndef ZM_ZM_DOWNDRAFT_PROPERTIES_IMPL_HPP
#define ZM_ZM_DOWNDRAFT_PROPERTIES_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_downdraft_properties. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_downdraft_properties(
  // Inputs
  const MemberType& team,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& jb, // updraft base level
  // Inputs/Outputs
  Int& jt, // updraft top level
  // Inputs
  const Int& j0, // level where updraft begins detraining
  // Inputs/Outputs
  Int& jd, // level of downdraft
  // Inputs
  const uview_1d<const Real>& z_int, // env altitude at interface
  const uview_1d<const Real>& dz, // layer thickness
  const uview_1d<const Real>& s_mid, // env dry energy of env [K] (normalized)
  const uview_1d<const Real>& q_mid, // env specific humidity
  const uview_1d<const Real>& h_env, // ambient env moist stat energy
  const uview_1d<const Real>& lambda, // fractional entrainment
  const Real& lambda_max, // fractional entrainment max
  const uview_1d<const Real>& qsthat, // interface interpolated qst
  const uview_1d<const Real>& hsthat, // interface interpolated hst
  const uview_1d<const Real>& gamhat, // interface interpolated gamma
  const uview_1d<const Real>& rprd, // rate of production of precip at that layer
  const uview_1d<const Real>& mflx_up, // updraft mass flux
  // Inputs/Outputs
  const uview_1d<Real>& mflx_dn, // downdraft mass flux
  const uview_1d<Real>& entr_dn, // downdraft entrainment rate
  const uview_1d<Real>& s_dnd, // dndraft dry energy [K] (normalized)
  const uview_1d<Real>& q_dnd, // dndraft specific humidity [kg/kg]
  const uview_1d<Real>& h_dnd, // dndraft moist energy
  const uview_1d<Real>& q_dnd_sat, // dndraft saturation specific humdity
  const uview_1d<Real>& evp, // evaporation rate
  Real& totevp) // total evap   for dndraft proportionality factor - see eq (4.106)
{
  //----------------------------------------------------------------------------
  // Purpose: Calculate properties of ZM downdrafts
  // Notes:
  // - Downward mass flux is scaled so that net flux (up-down) at cloud base in not negative
  // - No downdrafts if jd>=jb
  //----------------------------------------------------------------------------
  // Local variables
  Real ratmjb = 0; // ?
  //----------------------------------------------------------------------------
  // calculate downdraft mass flux
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    jt = ekat::impl::min(jt, jb-1);
    jd = ekat::impl::max(j0, jt+1);
    jd = ekat::impl::min(jd, jb);
    h_dnd(jd) = h_env(jd-1);
    if (jd < jb && lambda_max > 0) {
      // NOTE - this nonsensical lambda_max/lambda_max factor
      // was retained to preserve BFB results during ZM refactoring
      mflx_dn(jd) = -runtime_opt.alfa * lambda_max / lambda_max;
    }
  });
  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&](const Int& k) {
    if ((k > jd && k <= jb) && lambda_max > 0) {
      const Real dz_tmp = z_int(jd) - z_int(k);
      mflx_dn(k) = -runtime_opt.alfa / (2*lambda_max) * (std::exp(2*lambda_max*dz_tmp) - 1) / dz_tmp;
    }
  });
  team.team_barrier();

  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    if (lambda_max > 0 && jd < jb) {
      ratmjb = ekat::impl::min(std::abs(mflx_up(jb) / mflx_dn(jb)), Real(1));
    }
  });
  team.team_barrier();

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&](const Int& k) {
    if ((k >= jt && k <= jb) && lambda_max > 0 && jd < jb) {
      mflx_dn(k) *= ratmjb;
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // calculate downdraft entrainment and MSE
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    for (Int k = msg; k < pver; ++k) {
      if (k >= jt && lambda_max > 0) {
        entr_dn(k-1) = (mflx_dn(k-1) - mflx_dn(k)) / dz(k-1);
        const Real mdt = ekat::impl::min(mflx_dn(k), -ZMC::small);
        h_dnd(k) = (mflx_dn(k-1)*h_dnd(k-1) - dz(k-1)*entr_dn(k-1)*h_env(k-1)) / mdt;
      }
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // calculate downdraft specific humidity
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg+1, pver), [&](const Int& k) {
    if ((k >= jd && k <= jb) && lambda_max > 0 && jd < jb) {
      q_dnd_sat(k) = qsthat(k) + gamhat(k)*(h_dnd(k) - hsthat(k)) / (PC::LatVap.value*(1 + gamhat(k)));
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // downdraft quantities at source level
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    q_dnd(jd) = q_dnd_sat(jd);
    s_dnd(jd) = (h_dnd(jd) - PC::LatVap.value*q_dnd(jd)) / PC::Cpair.value;
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // calculate downdraft evaporation
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    for (Int k = msg+1; k < pver; ++k) {
      if (k >= jd && k < jb && lambda_max > 0) {
        q_dnd(k+1) = q_dnd_sat(k+1);
        evp(k) = -entr_dn(k)*q_mid(k) + (mflx_dn(k)*q_dnd(k) - mflx_dn(k+1)*q_dnd(k+1)) / dz(k);
        evp(k) = ekat::impl::max(evp(k), Real(0));
        const Real mdt = ekat::impl::min(mflx_dn(k+1), -ZMC::small);
        s_dnd(k+1) = ((PC::LatVap.value/PC::Cpair.value*evp(k) - entr_dn(k)*s_mid(k))*dz(k) + mflx_dn(k)*s_dnd(k)) / mdt;
        totevp -= dz(k)*entr_dn(k)*q_mid(k);
      }
    }
  });
  team.team_barrier();

  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    totevp += mflx_dn(jd)*q_dnd(jd) - mflx_dn(jb)*q_dnd(jb);
  });
}

} // namespace zm
} // namespace scream

#endif
