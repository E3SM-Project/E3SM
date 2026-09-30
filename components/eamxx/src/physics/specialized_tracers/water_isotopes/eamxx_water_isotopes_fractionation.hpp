#ifndef EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP
#define EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP

#include "share/core/eamxx_types.hpp"
#include "eamxx_water_isotopes_parameters.hpp"  // WaterIsotopologues, coefficient tables

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::exp/ekat::pow overloads for ekat::Pack
#include <ekat_kernel_assert.hpp>  // EKAT_KERNEL_REQUIRE_MSG

#include <cmath>  // std::exp/std::pow for the plain-Real path

namespace scream {
namespace wiso {

/*
 * Equilibrium isotopic fractionation factors for water isotopologues.
 *
 * Derived from the iCAM Fortran module water_isotopes.F90 (functions wiso_alpl
 * and wiso_alpi, original author David Noone). These are pure, device-callable
 * functions of temperature only; species is selected by a scalar enum (uniform
 * across a Pack, so no per-lane masking is needed for that). Templated on
 * ScalarT, which must be an ekat::Pack<Real,N> for some N -- N=1 stands in for
 * the scalar case, per EKAT's own convention (see e.g. p3_functions.hpp). exp()
 * and pow() are called unqualified so ADL selects the ekat::Pack overloads,
 * which exist for any N including 1.
 *
 * Convention: alpha *usually* constructed such that it is >=1, and 
 * thus follows:
 *   alpha = R_condensed / R_vapor  (>= 1),
 * i.e. the heavy isotope is preferentially retained in the condensed phase.
 * However, there are cases in the literature where alpha is defined differently,
 * and to avoid ambiguity adding a direction argument so that the source and
 * destination phases are always declared
 *
 * Both phases and elements share a single polynomial evaluator.
 */

// Which R-ratio the returned factor represents.
enum WisoAlphaDir {
  CondensedOverVapor = 0,  // raw table value R_condensed/R_vapor (>= 1)
  VaporOverCondensed = 1   // reciprocal, R_vapor/R_condensed (<= 1)
};

struct WaterIsotopeFractionation
{
private:
  // Mass-dependent scaling exponents
  static constexpr double H217O_exponent = 0.529;  // Schoenemann et al. (2014)
  static constexpr double HTO_exponent = 2.0;      // isoCAM3/5/6 assumption

  // Specify an arbitrarily low, non-Earth-system temperature to 
  // catch uninitialized memory or corrupted field.
  static constexpr double T_implausible = 50.0;  // [K]

  // In-range temperature substituted into dead lanes before the 1/T terms are
  // formed. 273.15K is valid for all alpha formulations.
  static constexpr double T_lane_fill = 273.15;  // [K]

  /* Temperature sanity, in two tiers.

     Tier 1 (always compiled in, including release): an implausible or NaN
     temperature ends the run.

     Tier 2 (debug builds only): a plausible temperature outside the fitted
     range means the polynomial is being extrapolated rather than evaluated.
     This may be the best we can do in some cases. Gated on NDEBUG because
     it fires per pack per level per step, which would swamp a production log.

     Lanes outside range_mask are ignored by both tiers: they hold padding, not
     data. */
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static void check_temperature(
      const ScalarT& t,
      const TemperatureBounds& b,
      const ekat::Mask<ScalarT::n>& range_mask,
      const char* caller)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    // isnan is tested separately: every comparison against NaN is false, so the
    // bounds test alone would let NaN through.
    const auto bad = (ekat::isnan(t) || t < RealT(T_implausible)) && range_mask;
    EKAT_KERNEL_REQUIRE_MSG(!bad.any(), caller);

#ifndef NDEBUG
    const auto extrapolating =
        (t < RealT(b.Tmin) || t > RealT(b.Tmax)) && range_mask;
    if (extrapolating.any()) {
      Kokkos::printf("WARNING: %s: T outside the fitted range [%g, %g] K;"
                     " extrapolating\n",
                     caller, double(b.Tmin), double(b.Tmax));
    }
#endif
  }

  // calculate 10^3 * ln(alpha) given a coefficient set and tempearture
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT ln_alpha_permil(const ScalarT& t, const PolynomialCoefficients& c)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    const ScalarT it  = RealT(1) / t;
    const ScalarT it2 = it * it;

    return ((RealT(c.T3)*t + RealT(c.T2))*t + RealT(c.T1))*t 
         + RealT(c.T0)
         + it*(RealT(c.T_1)
         + it*(RealT(c.T_2)
         + it*(RealT(c.T_3)
         + it*(RealT(c.T_4) + it2*RealT(c.T_6)))));
  }

  // Common fractionation logic given a species, temperature, direction
  template <typename ScalarT, typename BaseFunc>
  KOKKOS_INLINE_FUNCTION
  static ScalarT compute_alpha(const ScalarT& t,
                                const WaterIsotopologues species,
                                const WisoAlphaDir dir,
                                BaseFunc base_alpha)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    ScalarT alpha(1);

    switch (species) {
      case WaterIsotopologues::HDO:
        alpha = base_alpha(t, WaterIsotopologues::HDO);
        break;
      case WaterIsotopologues::H218O:
        alpha = base_alpha(t, WaterIsotopologues::H218O);
        break;
      case WaterIsotopologues::H217O:
        // Derived from H218O via mass-dependent fractionation
        alpha = pow(base_alpha(t, WaterIsotopologues::H218O), RealT(H217O_exponent));
        break;
      case WaterIsotopologues::HTO:
        // Derived from HDO via mass-dependent fractionation
        alpha = pow(base_alpha(t, WaterIsotopologues::HDO), RealT(HTO_exponent));
        break;
      case WaterIsotopologues::H216O:
      default:
        // Non-fractionating (alpha = 1)
        break;
    }

    // Apply direction
    return (dir == VaporOverCondensed) ? (RealT(1) / alpha) : alpha;
  }

public:
  // -----------------------------------------------------------------------
  // Equilibrium fractionation factor for either condensed phase.
  //
  //   alpha = exp( 1e-3 * (10^3 * ln alpha) )
  //
  // t is in Kelvin. The coefficient row is chosen by (phase, substituted
  // element) from the formulation held in `constants`; the 1e-3 undoes the
  // per-mil convention the source publications tabulate in.
  //
  // range_mask selects the lanes holding real data. Lanes outside it are
  // skipped by the temperature guard and replaced with an in-range value before
  // the 1/T terms are formed, so NaN padding left behind by upstream physics
  // cannot raise a spurious FPE. Callers holding a fully-populated pack (N=1
  // scalar callers always are) may use the overload that omits it.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_equilibrium(
      const ScalarT& t,
      const WaterIsotopologues species,
      const CondensedPhase phase,
      const WisoAlphaDir dir,
      const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params,
      const ekat::Mask<ScalarT::n>& range_mask = ekat::Mask<ScalarT::n>(true),
      const char* caller = "wiso::alpha_equilibrium")
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    using Params = WaterIsotopeParameters<RealT>;

    auto base = [&](const ScalarT& temp, WaterIsotopologues sp) -> ScalarT {
      const IsoElement el = Params::element_of(sp);

      // Bounds are per (phase, element): the default ice formulation draws its
      // two rows from two different studies with different fitted ranges.
      ScalarT t_live{RealT(T_lane_fill)};
      t_live.set(range_mask, temp);
      check_temperature(t_live, iso_params.tbounds(phase, el), range_mask, caller);

      return exp(RealT(1e-3) *
                 ln_alpha_permil(t_live, iso_params.alpha_eq_coeffs(phase, el)));
    };

    return compute_alpha(t, species, dir, base);
  }

  // -----------------------------------------------------------------------
  // Phase-specific spellings. These exist so call sites read as physics rather
  // than as table lookups; they add no logic.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_liquid_vapor(
      const ScalarT& t,
      const WaterIsotopologues species,
      const WisoAlphaDir dir,
      const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params,
      const ekat::Mask<ScalarT::n>& range_mask = ekat::Mask<ScalarT::n>(true))
  {
    return alpha_equilibrium(t, species, CondensedPhase::Liquid, dir, iso_params,
                             range_mask, "wiso::alpha_liquid_vapor");
  }

  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_ice_vapor(
      const ScalarT& t,
      const WaterIsotopologues species,
      const WisoAlphaDir dir,
      const WaterIsotopeParameters<typename ekat::ScalarTraits<ScalarT>::scalar_type>& iso_params,
      const ekat::Mask<ScalarT::n>& range_mask = ekat::Mask<ScalarT::n>(true))
  {
    return alpha_equilibrium(t, species, CondensedPhase::Ice, dir, iso_params,
                             range_mask, "wiso::alpha_ice_vapor");
  }
};

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP
