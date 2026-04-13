#ifndef ZM_ZM_CLOSURE_IMPL_HPP
#define ZM_ZM_CLOSURE_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_closure. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_closure(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Real& cape_threshold_in, // CAPE threshold for "cloud work function" (i.e. A)
  const Int& lcl, // index of lcl
  const Int& lel, // index of launch leve
  const Int& jt, // top of updraft
  const Int& mx, // base of updraft
  const Real& dsubcld, // thickness of subcloud layer
  const uview_1d<const Real>& z_mid, // altitude (m)
  const uview_1d<const Real>& z_int, // height of interface levels
  const uview_1d<const Real>& p_mid, // ambient pressure (mb)
  const uview_1d<const Real>& p_del, // pressure thickness of layers
  const uview_1d<const Real>& t_mid, // ambient temperature
  const uview_1d<const Real>& s_mid, // ambient dry static energy (normalized)
  const uview_1d<const Real>& q_mid, // ambient specific humidity
  const uview_1d<const Real>& qs, // ambient saturation specific humidity
  const uview_1d<const Real>& ql, // ambient liquid water mixing ratio
  const uview_1d<const Real>& s_int, // env. normalized dry static energy at intrfcs
  const uview_1d<const Real>& q_int, // environment specific humidity at interfaces
  const Real& t_pcl_lcl, // parcel temperature at LCL
  const uview_1d<const Real>& t_pcl, // parcel temperature
  const uview_1d<const Real>& q_pcl_sat, // parcel specific humidity
  const uview_1d<const Real>& s_upd, // updraft dry static energy (normalized)
  const uview_1d<const Real>& q_upd, // updraft specific humidity
  const uview_1d<const Real>& mflx_net, // net convective mass flux
  const uview_1d<const Real>& detr_up, // detrainment from updraft
  const uview_1d<const Real>& mflx_up, // updraft mass flux
  const uview_1d<const Real>& mflx_dn, // dndraft mass flux
  const uview_1d<const Real>& q_dnd, // dndraft specific humidity
  const uview_1d<const Real>& s_dnd, // dndraft dry static energy
  const Real& cape, // convective available potential energy
  // Outputs
  Real& cld_base_mass_flux) // cloud base mass flux
{
  //----------------------------------------------------------------------------
  // Purpose: calculate closure condition for ZM convection scheme using the
  //          revised quasi-equilibrium hypothesis of Z02, in which a
  //          quasi-equilibrium exists between the convective and large-scale
  //          modifications of the free-tropospheric CAPE, such that the net
  //          contribution is negilgible. This differs notably from AS74, where
  //          they assumed that CAPE changes from free-tropospheric and boundary
  //          layer changes are in balance. The Z02 revised closure is based on
  //          the observation that the total CAPE change is comparable to the
  //          CAPE change due to boundary layer thermodynamic changes.
  //----------------------------------------------------------------------------

  constexpr Real latvap = PC::LatVap.value;
  constexpr Real cpair  = PC::Cpair.value;
  constexpr Real grav   = PC::gravit.value;
  constexpr Real Rair   = PC::Rair.value;
  constexpr Real epsilo = PC::ep_2.value;

  // proportion to use liquid water from layer below
  constexpr Real beta = 0;

  // dtmdt: free tropospheric tendencies
  // dqmdt: free tropospheric tendencies
  // dboydt: integrand of cape change
  uview_1d<Real> dtmdt, dqmdt, dboydt;
  workspace.template take_many_contiguous_unsafe<3>(
    {"dtmdt", "dqmdt", "dboydt"},
    {&dtmdt, &dqmdt, &dboydt});

  Real dtbdt = 0; // sub-cloud layer tendencies
  Real dqbdt = 0; // sub-cloud layer tendencies
  Real dtldt = 0; // sub-cloud layer tendencies
  Real dadt  = 0; // CAPE consumption rate per unit cloud base mass flux (i.e. "F")

  //----------------------------------------------------------------------------
  // initialization
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&](const Int& k) {
    dtmdt(k)  = 0;
    dqmdt(k)  = 0;
    dboydt(k) = 0;
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Calculate sub-cloud tendencies of virtual temperature and humidity
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    cld_base_mass_flux = 0;
    Real eb   = p_mid(mx)*q_mid(mx) / (epsilo + q_mid(mx));
    dtbdt     = (1/dsubcld)
                * (mflx_up(mx)*(s_int(mx) - s_upd(mx))
                   + mflx_dn(mx)*(s_int(mx) - s_dnd(mx)));
    dqbdt     = (1/dsubcld)
                * (mflx_up(mx)*(q_int(mx) - q_upd(mx))
                   + mflx_dn(mx)*(q_int(mx) - q_dnd(mx)));
    Real debdt = epsilo*p_mid(mx) / ((epsilo + q_mid(mx))*(epsilo + q_mid(mx))) * dqbdt;
    Real log_arg = 3.5*std::log(t_mid(mx)) - std::log(eb) - 4.805;
    dtldt     = -2840 * (3.5/t_mid(mx)*dtbdt - debdt/eb) / (log_arg*log_arg);
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Calculate dtmdt & dqmdt
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver-1), [&](const Int& k) {
    // cloud top
    if (k == jt) {
      dtmdt(k) = (1/p_del(k))
                 * (mflx_up(k+1)*(s_upd(k+1) - s_int(k+1) - latvap/cpair*ql(k+1))
                    + mflx_dn(k+1)*(s_dnd(k+1) - s_int(k+1)));
      dqmdt(k) = (1/p_del(k))
                 * (mflx_up(k+1)*(q_upd(k+1) - q_int(k+1) + ql(k+1))
                    + mflx_dn(k+1)*(q_dnd(k+1) - q_int(k+1)));
    }
    // below cloud top
    if (k > jt && k < mx) {
      dtmdt(k) = (mflx_net(k  )*(s_int(k  ) - s_mid(k))
                  + mflx_net(k+1)*(s_mid(k  ) - s_int(k+1))) / p_del(k)
                 - latvap/cpair * detr_up(k)*(beta*ql(k) + (1-beta)*ql(k+1));
      dqmdt(k) = (mflx_up(k+1)*(q_upd(k+1) - q_int(k+1) + cpair/latvap*(s_upd(k+1) - s_mid(k)))
                  - mflx_up(k  )*(q_upd(k  ) - q_int(k  ) + cpair/latvap*(s_upd(k  ) - s_mid(k)))
                  + mflx_dn(k+1)*(q_dnd(k+1) - q_int(k+1) + cpair/latvap*(s_dnd(k+1) - s_mid(k)))
                  - mflx_dn(k  )*(q_dnd(k  ) - q_int(k  ) + cpair/latvap*(s_dnd(k  ) - s_mid(k)))) / p_del(k)
                 + detr_up(k)*(beta*ql(k) + (1-beta)*ql(k+1));
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Calculate dboydt (integrand of cape change)
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, msg, pver), [&](const Int& k) {
    // levels between parcel launch and LCL
    if (k > lcl && k < mx) {
      Real thetavp = t_pcl(k) * std::pow(Real(1000)/p_mid(k), Rair/cpair)
                     * (1 + 0.608*q_mid(mx));
      Real thetavm = t_mid(k) * std::pow(Real(1000)/p_mid(k), Rair/cpair)
                     * (1 + 0.608*q_mid(k));
      dboydt(k) = (dtbdt/t_mid(mx) + 0.608/(1 + 0.608*q_mid(mx))*dqbdt
                   - dtmdt(k)/t_mid(k) - 0.608/(1 + 0.608*q_mid(k))*dqmdt(k))
                  * grav * thetavp/thetavm;
    }
    // levels between LCL and cloud top
    if (k >= lel && k <= lcl) {
      Real thetavp = t_pcl(k) * std::pow(Real(1000)/p_mid(k), Rair/cpair)
                     * (1 + 1.608*q_pcl_sat(k) - q_mid(mx));
      Real thetavm = t_mid(k) * std::pow(Real(1000)/p_mid(k), Rair/cpair)
                     * (1 + 0.608*q_mid(k));
      Real dqsdtp  = q_pcl_sat(k) * (1 + q_pcl_sat(k)/epsilo) * epsilo*latvap
                     / (Rair*t_pcl(k)*t_pcl(k));
      Real dtpdt   = t_pcl(k) / (1 + latvap/cpair*(dqsdtp - q_pcl_sat(k)/t_pcl(k)))
                     * (dtbdt/t_mid(mx) + latvap/cpair*(dqbdt/t_pcl_lcl
                        - q_mid(mx)/(t_pcl_lcl*t_pcl_lcl)*dtldt));
      dboydt(k) = ((dtpdt/t_pcl(k) + 1/(1 + 1.608*q_pcl_sat(k) - q_mid(mx))
                    * (1.608*dqsdtp*dtpdt - dqbdt))
                   - (dtmdt(k)/t_mid(k) + 0.608/(1 + 0.608*q_mid(k))*dqmdt(k)))
                  * grav * thetavp/thetavm;
    }
  });
  team.team_barrier();

  //----------------------------------------------------------------------------
  // vertically integrate buoyancy change
  // kmin = lel, kmax = mx - 1 (TeamVectorRange is half-open: [lel, mx))
  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, lel, mx),
    [&](const Int& k, Real& dadt_sum) {
      dadt_sum += dboydt(k) * (z_int(k) - z_int(k+1));
    }, dadt);
  team.team_barrier();

  //----------------------------------------------------------------------------
  // Calculate cloud base mass flux - see eq (8) in Z02
  Kokkos::single(Kokkos::PerTeam(team), [&]() {
    Real dltaa = -1 * (cape - cape_threshold_in);
    if (dadt != 0) {
      cld_base_mass_flux = Kokkos::max(dltaa/runtime_opt.tau/dadt, Real(0));
    }
  });

  workspace.template release_many_contiguous<3>({&dtmdt, &dqmdt, &dboydt});
}

} // namespace zm
} // namespace scream

#endif
