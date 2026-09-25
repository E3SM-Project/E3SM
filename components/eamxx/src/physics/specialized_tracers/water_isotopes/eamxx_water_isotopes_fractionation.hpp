#ifndef EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP
#define EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP

#include "share/core/eamxx_types.hpp"  // for scream::Real and scream::sp()
#include "eamxx_water_isotopes_constants.hpp"  // WaterIsotopologues, coefficient tables

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::exp/ekat::pow overloads for ekat::Pack (found via ADL)
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
 * ScalarT so they work for both a plain Real and an ekat::Pack<Real,N>: exp()
 * and pow() are called unqualified so ADL selects the ekat::Pack overloads for
 * packs and std for plain scalars (matching PhysicsFunctions::exner_function).
 *
 * Convention: alpha *usually* constructed such that it is >=1, and 
 * the follow this convention:
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

namespace impl {

/* Pairs ScalarT with its "which lanes are live" type: ekat::Mask<N> for a
   Pack<T,N>, a plain bool for a bare scalar. This lets the temperature guard
   and the padded-lane neutralization below be written once and used from both
   the packed physics path and the scalar unit tests. */
template <typename ScalarT>
struct LaneTraits {
  using mask_type = bool;

  KOKKOS_INLINE_FUNCTION static mask_type all_lanes() { return true; }
  KOKKOS_INLINE_FUNCTION static bool any(const mask_type& m) { return m; }
  KOKKOS_INLINE_FUNCTION static mask_type is_nan(const ScalarT& v) {
    return ekat::impl::is_nan(v);
  }
  KOKKOS_INLINE_FUNCTION
  static ScalarT select(const ScalarT& v, const mask_type& live, const ScalarT& fill) {
    return live ? v : fill;
  }
};

template <typename T, int N>
struct LaneTraits<ekat::Pack<T,N>> {
  using mask_type = ekat::Mask<N>;

  KOKKOS_INLINE_FUNCTION static mask_type all_lanes() { return mask_type(true); }
  KOKKOS_INLINE_FUNCTION static bool any(const mask_type& m) { return m.any(); }
  KOKKOS_INLINE_FUNCTION static mask_type is_nan(const ekat::Pack<T,N>& v) {
    return ekat::isnan(v);
  }
  KOKKOS_INLINE_FUNCTION
  static ekat::Pack<T,N> select(const ekat::Pack<T,N>& v, const mask_type& live,
                                const ekat::Pack<T,N>& fill) {
    ekat::Pack<T,N> out(fill);
    out.set(live, v);
    return out;
  }
};

} // namespace impl

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
      const typename impl::LaneTraits<ScalarT>::mask_type& range_mask,
      const char* caller)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    using LT    = impl::LaneTraits<ScalarT>;

    // isnan is tested separately: every comparison against NaN is false, so the
    // bounds test alone would let NaN through.
    const auto bad = (LT::is_nan(t) || t < RealT(T_implausible)) && range_mask;
    EKAT_KERNEL_REQUIRE_MSG(!LT::any(bad), caller);

#ifndef NDEBUG
    const auto extrapolating =
        (t < RealT(b.Tmin) || t > RealT(b.Tmax)) && range_mask;
    if (LT::any(extrapolating)) {
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
  // cannot raise a spurious FPE. Callers holding a fully-populated pack (or a
  // bare scalar) may use the overload that omits it.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_equilibrium(
      const ScalarT& t,
      const WaterIsotopologues species,
      const CondensedPhase phase,
      const WisoAlphaDir dir,
      const WaterIsotopeConstants<typename ekat::ScalarTraits<ScalarT>::scalar_type>& constants,
      const typename impl::LaneTraits<ScalarT>::mask_type& range_mask =
          impl::LaneTraits<ScalarT>::all_lanes(),
      const char* caller = "wiso::alpha_equilibrium")
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    using Constants = WaterIsotopeConstants<RealT>;
    using LT = impl::LaneTraits<ScalarT>;

    auto base = [&](const ScalarT& temp, WaterIsotopologues sp) -> ScalarT {
      const IsoElement el = Constants::element_of(sp);

      // Bounds are per (phase, element): the default ice formulation draws its
      // two rows from two different studies with different fitted ranges.
      check_temperature(temp, constants.tbounds(phase, el), range_mask, caller);

      const ScalarT t_live = LT::select(temp, range_mask, ScalarT(RealT(T_lane_fill)));

      return exp(RealT(1e-3) *
                 ln_alpha_permil(t_live, constants.alpha_eq_coeffs(phase, el)));
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
      const WaterIsotopeConstants<typename ekat::ScalarTraits<ScalarT>::scalar_type>& constants,
      const typename impl::LaneTraits<ScalarT>::mask_type& range_mask =
          impl::LaneTraits<ScalarT>::all_lanes())
  {
    return alpha_equilibrium(t, species, CondensedPhase::Liquid, dir, constants,
                             range_mask, "wiso::alpha_liquid_vapor");
  }

  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_ice_vapor(
      const ScalarT& t,
      const WaterIsotopologues species,
      const WisoAlphaDir dir,
      const WaterIsotopeConstants<typename ekat::ScalarTraits<ScalarT>::scalar_type>& constants,
      const typename impl::LaneTraits<ScalarT>::mask_type& range_mask =
          impl::LaneTraits<ScalarT>::all_lanes())
  {
    return alpha_equilibrium(t, species, CondensedPhase::Ice, dir, constants,
                             range_mask, "wiso::alpha_ice_vapor");
  }
};

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP
