#ifndef P3_EVAPORATE_RAIN_IMPL_HPP
#define P3_EVAPORATE_RAIN_IMPL_HPP

#include "p3_functions.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_evap_tscale_weight(const Spack& dt_over_tau, Spack& weight, const Smask& context)
{
  /*
    Returns weighting between 0 and 1 for how much of the instantaneous
    evaporation rate and how much of the equilibrium evaporation rate to
    blend to get the timestep-average rain evaporation rate
  */

  //weight.set(context, (1 - exp(-dt_over_tau) )/dt_over_tau );
  weight.set(context, -expm1(-dt_over_tau)/dt_over_tau );

} //end tscale_weight

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_evap_equilib_tend(const Spack& A_c,const Spack& ab,const Spack& tau_eff,
			 const Spack& tau_r, Spack& tend, const Smask& context)
{
  /*
    In equilibrium, the total evaporation must balance the tendency A_c from
    all other processes. The rain evaporation is the fraction (1/tau_r)/(1/tau_eff)
    of the total tendency and ab corrects for saturation changes due to evaporative
    cooling.
  */

  //Sign convention: Negative A_c causes a supersaturation deficit which needs to be removed
  //by evaporation (which is signed positive when active) to maintain equilibrium. Thus
  //A_c needs a negative sign here since other terms are always positive.
  tend.set(context,-A_c/ab*tau_eff/tau_r);

} //end equilib_tend

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_evap_instant_tend(const Spack& ssat_r, const Spack& ab, const Spack& tau_r,
			 Spack& tend, const Smask& context)
{
  /*
    The instantaneous rain evap tendency is just the absolute supersaturation
    ssat_r divided by the supersaturation removal timescale for rain tau_r
    corrected for the effect of evaporative cooling on saturation ab.
  */

  //sign convention: ssat_r must be <0 for evap, other terms are always positive,
  //and we want evap rate positive... so put a minus sign in front.
  tend.set(context,-ssat_r/(ab*tau_r) );

}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::evaporate_rain(
  const Spack& qr_incld, const Spack& qc_incld, const Spack& nr_incld, const Spack& qi_incld,
  const Spack& cld_frac_l, const Spack& cld_frac_r, const Spack& qv, const Spack& qv_prev,
  const Spack& qv_sat_l, const Spack& qv_sat_i, const Spack& ab, const Spack& abi,
  const Spack& epsr, const Spack& epsi_tot, const Spack& t_atm, const Spack& t_atm_prev,
  const Spack& latent_heat_sublim, const Spack& dqsdt, const Scalar& dt,
  Spack& qr2qv_evap_tend, Spack& nr_evap_tend,
  const Smask& context)
{
  /* Evaporation is basically (qv - sv_sat)/(tau_eff*ab) where tau_eff
     is the total effective supersaturation removal timescale
     and ab is the psychrometric correction for condensational heating
     changing qv_sat. This formulation depends sensitively on ssat_r, which
     can change rapidly within a timestep because liquid saturation
     adjustment has a relaxation timescale of seconds. For accuracy and
     stability, we analytically integrate ssat_r over the timestep under
     the simplifying assumption that all processes other than saturation
     relaxation are a constant source/sink term A_c. See Morrison+Milbrandt 2015
     https://doi.org/10.1175/JAS-D-14-0065.1 and Morrison+Grabowski 2008
     https://doi.org/10.1175/2007JAS2374.1 for details. */

  //Initialize variables
  qr2qv_evap_tend = 0;
  nr_evap_tend = 0;
  const Scalar inv_dt = 1/dt;
  constexpr Scalar QSMALL   = C::QSMALL;
  constexpr Scalar Tmelt  = C::Tmelt;
  constexpr Scalar inv_cp = 1/C::Cpair;

  //Compute absolute supersaturation.
  //Ignore the difference between clear-sky and cell-ave qv and T
  //because micro lacks the info to reliably reconstruct macrophys
  //subgrid variability
  Spack ssat_r = qv - qv_sat_l;

  //Cloud fraction in clear-sky conditions has been set to mincld
  //to avoid divide-by-zero problems. Because rain evap only happens
  //in rainy portions outside cloud, setting clear-sky cloud fraction
  //to mincld reduces evaporating area. We fix that here by computing
  //a temporary cloud fraction which is zero if cloud condensate is small.
  Spack cld_frac;
  const auto set_cld_frac_zero = (qc_incld + qi_incld) < sp(1.e-6) && context;
  cld_frac.set(set_cld_frac_zero, 0);
  cld_frac.set(!set_cld_frac_zero && context, cld_frac_l);

  //Only evaporate in the rainy area outside cloud when subsaturated
  //Note: ignoring case where cell initially supersaturated but other
  //processes would make it subsaturated within 1 timestep.
  const Smask qr_ge_qsmall = qr_incld >= QSMALL;
  const Smask is_subsat = ssat_r < 0;
  const Smask is_evap_area = cld_frac_r > cld_frac;
  const Smask is_rain_evap = qr_ge_qsmall && is_subsat && is_evap_area && context;
  if (is_rain_evap.any()){

    //if qr_incld<QSMALL, epsr=0 causes div by 0 error for tau_r even though it isn't used.
    Spack tau_r(0);
    tau_r.set(is_rain_evap, 1/epsr);

    //Compute total effective inverse saturation removal timescale eps_eff
    //qc saturation is handled by macrophysics so the qc saturation removal timescale is
    //not included here. Below freezing, eps_eff is the sum of the inverse saturation
    //removal timescales for liquid and ice. The ice term has extra scaling terms to convert
    //it from being relative to ice to liquid instead. Also compute the constant source/sink
    //term A_c for analytic integration. See Eq C3 and C4 of Morrison+Milbrandt 2015
    //https://doi.org/10.1175/JAS-D-14-0065.1, respectively.
    const Smask is_freezing = t_atm < Tmelt && context;
    const Smask not_freezing = !is_freezing && context;
    Spack eps_eff, A_c;
    if (is_freezing.any()){
      eps_eff.set(is_freezing,epsr + epsi_tot*(1 + latent_heat_sublim*inv_cp*dqsdt)/abi);
      A_c.set(is_freezing,(qv - qv_prev)*inv_dt - dqsdt*(t_atm-t_atm_prev)*inv_dt
	      - (qv_sat_l - qv_sat_i)*(1 + latent_heat_sublim*inv_cp*dqsdt)/abi*epsi_tot );
    }
    if (not_freezing.any()){
      eps_eff.set(not_freezing,epsr);
      A_c.set(not_freezing, (qv - qv_prev)*inv_dt - dqsdt*(t_atm-t_atm_prev)*inv_dt );
    }

    //Set lower bound on eps_eff to prevent division by zero
    eps_eff.set(eps_eff<1e-20 && context, 1e-20);
    const Spack tau_eff = 1/eps_eff;

    //If qr is posive but tiny, evap all qr if subsaturated at all.
    Smask is_qr_tiny = qr_incld < 1e-12 && qv/qv_sat_l < 0.999;
    const Smask not_qr_tiny = !is_qr_tiny && is_rain_evap;
    is_qr_tiny = is_qr_tiny && is_rain_evap;
    if (is_qr_tiny.any()){
      qr2qv_evap_tend.set(is_qr_tiny,qr_incld*inv_dt );
    }

    //If qr is reasonably big, compute timestep-averaged evap rate as the weighted ave of
    //instantaneous and equilibrium evap rates with weighting timescale tscale_weight. L'Hospital's
    //rull shows tscale_weight is 1 in the limit of small dt. It approaches 0 as dt gets big.
    if (not_qr_tiny.any()){
      Spack tscale_weight, equilib_tend, instant_tend;
      rain_evap_tscale_weight(dt/tau_eff,tscale_weight,is_rain_evap);
      rain_evap_equilib_tend(A_c,ab,tau_eff,tau_r,equilib_tend,is_rain_evap);
      rain_evap_instant_tend(ssat_r, ab, tau_r,instant_tend,is_rain_evap);

      qr2qv_evap_tend.set(not_qr_tiny,
			  instant_tend*tscale_weight
			  + equilib_tend*(1-tscale_weight) );

    }

    //Limit evap from exceeding saturation deficit. Analytic integration
    //would prevent this from happening if A_c was part of microphysics
    //timestepping, but it isn't.
    const Smask is_overevap=qr2qv_evap_tend > -ssat_r*inv_dt/ab && is_rain_evap;
    qr2qv_evap_tend.set(is_overevap, -ssat_r*inv_dt/ab );

    //To maintain equilibrium, the equilibrium evaporation tendency must be
    //negative (adding mass) if A_c (other processes) are losing mass. We don't
    //allow rain evap to also condense by forcing qr2qv_evap_tend to be positive
    qr2qv_evap_tend.set(is_rain_evap && (qr2qv_evap_tend<0), 0);

    //We can't evaporate more rain mass than we had to start with
    //Note: We're applying to rainy region outside cloud here because
    //qr inside cloud should be protected from evap. Conversion to rainy-area
    //average just below scales by (cldfrac_r - cld_frac)/cldfrac_r < 1 so
    //total qr isn't pushed negative.
    qr2qv_evap_tend.set(is_rain_evap && (qr2qv_evap_tend>qr_incld*inv_dt),
			qr_incld*inv_dt);

    //Evap rate so far is an average over the rainy area outside clouds.
    //Turn this into an average over the entire raining area
    qr2qv_evap_tend.set(is_rain_evap, qr2qv_evap_tend*(cld_frac_r-cld_frac)/cld_frac_r);

    //Let nr remove drops proportionally to mass change
    nr_evap_tend.set(is_rain_evap, qr2qv_evap_tend*(nr_incld/qr_incld));

  } //end if (rain_evap.any()
} //end evaporate rain


} // namespace p3
} // namespace scream

#endif
