#ifndef ZM_ZM_CONV_EVAP_IMPL_HPP
#define ZM_ZM_CONV_EVAP_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_conv_evap. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_conv_evap(
  // Inputs
  const MemberType& team,
  const ZmRuntimeOpt& runtime_opt,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Real& time_step, // model time step                         [s]
  const uview_1d<const Real>& p_mid, // midpoint pressure                       [Pa]
  const uview_1d<const Real>& p_del, // layer thickness                         [Pa]
  const uview_1d<const Real>& t_mid, // temperature                             [K]
  const uview_1d<const Real>& q_mid, // water vapor                             [kg/kg]
  const uview_1d<const Real>& prdprec, // precipitation production                [kg/kg/s]
  const uview_1d<const Real>& cldfrc, // cloud fraction
  // Inputs/Outputs
  const uview_1d<Real>& tend_s, // heating rate                            [J/kg/s]
  const uview_1d<Real>& tend_q, // water vapor tendency                    [kg/kg/s]
  // Outputs
  const uview_1d<Real>& tend_s_snwprd, // Heating rate of snow production         [J/kg/s]
  const uview_1d<Real>& tend_s_snwevmlt, // Heating rate of snow evap/melt          [J/kg/s]
  // Inputs/Outputs
  Real& prec, // Convective-scale prec rate              [m/s]
  // Outputs
  Real& snow, // Convective-scale snow rate              [m/s]
  const uview_1d<Real>& ntprprd, // net precip production in layer          [kg/kg/s]
  const uview_1d<Real>& ntsnprd, // net snow production in layer            [kg/kg/s]
  const uview_1d<Real>& flxprec, // Convective flux of prec at interfaces   [kg/m2/s]
  const uview_1d<Real>& flxsnow) // Convective flux of snow at interfaces   [kg/m2/s])
{
  // perturbation constant used in pergro snow fraction guard
  const Real pergro_perturbation = 8.64e-11;

#ifdef PERGRO
  const bool pergro_active = true;
#else
  const bool pergro_active = false;
#endif

  // cldfrc_fice temperature thresholds for convective snow fraction
  // (inlined from cloud_fraction.F90 cldfrc_fice, snow-fraction branch)
  const Real tmax_fsnow = PC::Tmelt.value;        // max temperature for transition to convective snow
  const Real tmin_fsnow = PC::Tmelt.value - 5;  // min temperature for transition to convective snow

  // The entire k-loop is sequential due to top-to-bottom flux dependencies
  // (flxprec and flxsnow at interface k+1 depend on interface k)
  Kokkos::single(Kokkos::PerTeam(team), [&] () {

    // convert input precip to kg/m2/s
    prec = prec * 1000;

    // determine saturation vapor pressure
    // (qs_k computed per level inside the loop via qsat_hPa)

    // determine ice fraction in rain production
    // (fsnow_conv_k computed per level inside the loop via inlined cldfrc_fice logic)

    // zero the flux integrals on the top boundary
    flxprec(0) = 0;
    flxsnow(0) = 0;
    Real evpvint = 0;

    for (int k = 0; k < pver; ++k) {

      // determine saturation vapor pressure and mixing ratio
      Real es_k, qs_k;
      qsat(t_mid(k), p_mid(k), runtime_opt, es_k, qs_k);

      // determine convective snow fraction for this level (from cldfrc_fice)
      // If warmer than tmax_fsnow then water phase
      Real fsnow_conv_k;
      if (t_mid(k) > tmax_fsnow) {
        fsnow_conv_k = 0;
      // If colder than tmin_fsnow then ice phase
      } else if (t_mid(k) < tmin_fsnow) {
        fsnow_conv_k = 1;
      // Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
      } else {
        fsnow_conv_k = (tmax_fsnow - t_mid(k)) / (tmax_fsnow - tmin_fsnow);
      }

      // snow production rate from microphysics; 0 when zm_microp not active
      const Real prdsnow_k = 0;

      // Melt snow falling into layer, if necessary.
      Real flxsntm, snowmlt;
      if (runtime_opt.old_snow) {
        if (t_mid(k) > PC::Tmelt.value) {
          flxsntm = 0;
          snowmlt = flxsnow(k) * PC::gravit.value / p_del(k);
        } else {
          flxsntm = flxsnow(k);
          snowmlt = 0;
        }
      } else {
        if (t_mid(k) > PC::Tmelt.value) {
          Real dum = -PC::LatIce.value / PC::Cpair.value * flxsnow(k) * PC::gravit.value / p_del(k) * time_step;
          if (t_mid(k) + dum <= PC::Tmelt.value) {
            dum = (t_mid(k) - PC::Tmelt.value) * PC::Cpair.value / PC::LatIce.value / time_step;
            dum = dum / (flxsnow(k) * PC::gravit.value / p_del(k));
            dum = ekat::impl::max(Real(0), dum);
            dum = ekat::impl::min(Real(1), dum);
          } else {
            dum = 1;
          }
          dum = dum * ZMC::omsm;
          flxsntm = flxsnow(k) * (1 - dum);
          snowmlt = dum * flxsnow(k) * PC::gravit.value / p_del(k);
        } else {
          flxsntm = flxsnow(k);
          snowmlt = 0;
        }
      }

      // relative humidity depression must be > 0 for evaporation
      Real evplimit = ekat::impl::max(1 - q_mid(k) / qs_k, Real(0));

      // total evaporation depends on flux in the top of the layer
      Real evpprec = runtime_opt.ke * (1 - cldfrc(k)) * evplimit * std::sqrt(flxprec(k));

      // Don't let evaporation supersaturate layer (approx).
      evplimit = ekat::impl::max(Real(0), (qs_k - q_mid(k)) / time_step);

      // Don't evaporate more than is falling into the layer from above.
      evplimit = ekat::impl::min(evplimit, flxprec(k) * PC::gravit.value / p_del(k));

      // Total evaporation cannot exceed input precipitation
      evplimit = ekat::impl::min(evplimit, (prec - evpvint) * PC::gravit.value / p_del(k));

      evpprec = ekat::impl::min(evplimit, evpprec);

      if (!runtime_opt.old_snow) {
        evpprec = ekat::impl::max(Real(0), evpprec);
        evpprec = evpprec * ZMC::omsm;
      }

      // evaporation of snow depends on snow fraction
      Real evpsnow;
      if (flxprec(k) > 0) {
        Real work1 = ekat::impl::min(ekat::impl::max(Real(0), flxsntm / flxprec(k)), Real(1));
        if (!runtime_opt.old_snow && prdsnow_k > prdprec(k)) work1 = 1;
        evpsnow = evpprec * work1;
      } else {
        evpsnow = 0;
      }

      // vertically integrated evaporation
      evpvint = evpvint + evpprec * p_del(k) / PC::gravit.value;

      // net precip production
      ntprprd(k) = prdprec(k) - evpprec;

      // net snow production
      if (runtime_opt.old_snow) {
        Real work1;
        if (pergro_active) {
          work1 = ekat::impl::min(ekat::impl::max(Real(0), flxsnow(k) / (flxprec(k) + pergro_perturbation)), Real(1));
        } else {
          if (flxprec(k) > 0) {
            work1 = ekat::impl::min(ekat::impl::max(Real(0), flxsnow(k) / flxprec(k)), Real(1));
          } else {
            work1 = 0;
          }
        }
        Real work2 = ekat::impl::max(fsnow_conv_k, work1);
        if (snowmlt > 0) work2 = 0;
        ntsnprd(k) = prdprec(k) * work2 - evpsnow - snowmlt;
        tend_s_snwprd  (k) = prdprec(k) * work2 * PC::LatIce.value;
        tend_s_snwevmlt(k) = -(evpsnow + snowmlt) * PC::LatIce.value;
      } else {
        ntsnprd(k) = prdsnow_k - ekat::impl::min(flxsnow(k) * PC::gravit.value / p_del(k), evpsnow + snowmlt);
        tend_s_snwprd  (k) = prdsnow_k * PC::LatIce.value;
        tend_s_snwevmlt(k) = -ekat::impl::min(flxsnow(k) * PC::gravit.value / p_del(k), evpsnow + snowmlt) * PC::LatIce.value;
      }

      // precipitation fluxes
      flxprec(k+1) = flxprec(k) + ntprprd(k) * p_del(k) / PC::gravit.value;
      flxsnow(k+1) = flxsnow(k) + ntsnprd(k) * p_del(k) / PC::gravit.value;

      flxprec(k+1) = ekat::impl::max(flxprec(k+1), Real(0));
      flxsnow(k+1) = ekat::impl::max(flxsnow(k+1), Real(0));

      // heating and moistening due to evaporation
      if (runtime_opt.old_snow) {
        tend_s(k) = -evpprec * PC::LatVap.value + ntsnprd(k) * PC::LatIce.value;
      } else {
        tend_s(k) = -evpprec * PC::LatVap.value + tend_s_snwevmlt(k);
      }
      tend_q(k) = evpprec;

    } // k loop

    // protect against rounding error
    if (!runtime_opt.old_snow) {
      if (flxsnow(pver) > flxprec(pver)) {
        Real dum = (flxsnow(pver) - flxprec(pver)) * PC::gravit.value;
        for (int k = pver - 1; k >= 0; --k) {
          if (ntsnprd(k) > ntprprd(k) && dum > 0) {
            ntsnprd(k) = ntsnprd(k) - dum / p_del(k);
            tend_s_snwevmlt(k) = tend_s_snwevmlt(k) - dum / p_del(k) * PC::LatIce.value;
            tend_s(k) = tend_s(k) - dum / p_del(k) * PC::LatIce.value;
            dum = 0;
          }
        }
        flxsnow(pver) = flxprec(pver);
      }
    }

    // set output precipitation rates (m/s)
    prec = flxprec(pver) / 1000;
    snow = flxsnow(pver) / 1000;

  }); // Kokkos::single
}

} // namespace zm
} // namespace scream

#endif
