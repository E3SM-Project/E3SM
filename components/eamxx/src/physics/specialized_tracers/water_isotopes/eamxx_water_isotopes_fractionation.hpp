#ifndef EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP
#define EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP

#include "share/core/eamxx_types.hpp"  // for scream::Real and scream::sp()
#include "eamxx_water_isotopes_constants.hpp"  // isotopic constants (available for future use)

#include <ekat_pack.hpp>
#include <ekat_pack_math.hpp>  // ekat::exp/ekat::pow overloads for ekat::Pack (found via ADL)

#include <cmath>  // std::exp/std::pow for the plain-Real path

namespace scream {
namespace wiso {

/*
 * Equilibrium isotopic fractionation factors for water isotopologues.
 *
 * Ported from the iCAM Fortran module water_isotopes.F90 (functions wiso_alpl
 * and wiso_alpi, author David Noone). These are pure, device-callable functions
 * of temperature only; species is selected by a scalar enum (uniform across a
 * Pack, so no per-lane masking is needed). Templated on ScalarT so they work
 * for both a plain Real and an ekat::Pack<Real,N>: exp() and pow() are called
 * unqualified so ADL selects the ekat::Pack overloads for packs and std for
 * plain scalars (matching PhysicsFunctions::exner_function).
 *
 * Convention: the underlying tables return the vapor->condensed enrichment
 *   alpha = R_condensed / R_vapor  (>= 1),
 * i.e. the heavy isotope is preferentially retained in the condensed phase.
 * The desired direction is chosen explicitly via WisoAlphaDir (required
 * argument) so every call site states its intent.
 *
 * Scope: EQUILIBRIUM factors only. Kinetic / supersaturation effects (iCAM
 * wiso_akel / wiso_akci, active below ~253 K) are intentionally NOT included
 * here; see the "KINETIC HOOK" note below for where a future
 * alpha_ice_vapor_kinetic(t, species, s_ice, ...) sibling would attach.
 *
 * Note: Isotopic constants are now centralized in eamxx_water_isotopes_constants.hpp
 * and available for other modules. The coefficients below are embedded for
 * backwards compatibility; future work may refactor to use the centralized constants.
 */

// Water isotopologues. HDO and H218O are computed directly from the tables;
// H217O and HTO are derived by mass-dependent power laws. H216O (ordinary
// water) is non-fractionating (alpha == 1).
enum WisoSpecies {
  H216O = 0,  // ordinary water; alpha == 1
  HDO   = 1,  // HD16O (deuterium)
  H218O = 2,  // H2-18O (oxygen-18)
  H217O = 3,  // H2-17O; = alpha(H218O)^0.529 (Schoenemann et al. 2014)
  HTO   = 4   // HT16O (tritiated water); = alpha(HDO)^2.0 (isoCAM3)
};

// Which R-ratio the returned factor represents.
enum WisoAlphaDir {
  CondensedOverVapor = 0,  // raw table value R_condensed/R_vapor (>= 1)
  VaporOverCondensed = 1   // reciprocal, R_vapor/R_condensed (<= 1)
};

struct WaterIsotopeFractionation
{
  using Real = scream::Real;

  // -----------------------------------------------------------------------
  // Liquid <-> vapor equilibrium fractionation factor.
  //
  // Horita & Wesolowski (1994). Two functional forms by species; T in Kelvin.
  //   HDO:   alpha = exp( a*T^3 + b*T^2 + c*T + d + e/T^3 )
  //   H218O: alpha = exp( a/T^3 + b/T^2 + c/T + d )
  // Coefficients transcribed directly from water_isotopes.F90 (lines 167-171).
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_liquid_vapor(const ScalarT& t,
                                    const WisoSpecies species,
                                    const WisoAlphaDir dir)
  {
    // HDO (Horita & Wesolowski 1994): a*T^3 + b*T^2 + c*T + d + e/T^3
    constexpr Real hdo_a = sp(1158.8e-12);
    constexpr Real hdo_b = sp(-1620.1e-9);
    constexpr Real hdo_c = sp(794.84e-6);
    constexpr Real hdo_d = sp(-161.04e-3);
    constexpr Real hdo_e = sp(2.9992e6);

    // H2-18O (Horita & Wesolowski 1994): a/T^3 + b/T^2 + c/T + d
    constexpr Real o18_a = sp(0.35041e6);
    constexpr Real o18_b = sp(-1.6664e3);
    constexpr Real o18_c = sp(6.7123);
    constexpr Real o18_d = sp(-7.685e-3);

    ScalarT alpha(1);

    switch (species) {
      case HDO: {
        const ScalarT t2 = t * t;
        const ScalarT t3 = t2 * t;
        alpha = exp(hdo_a * t3 + hdo_b * t2 + hdo_c * t + hdo_d + hdo_e / t3);
        break;
      }
      case H218O: {
        const ScalarT it  = sp(1) / t;
        const ScalarT it2 = it * it;
        const ScalarT it3 = it2 * it;
        alpha = exp(o18_a * it3 + o18_b * it2 + o18_c * it + o18_d);
        break;
      }
      case H217O:
        // Mass-dependent from H2-18O (Schoenemann et al. 2014).
        alpha = pow(alpha_liquid_vapor(t, H218O, CondensedOverVapor), sp(0.529));
        break;
      case HTO:
        // From HDO (isoCAM3).
        alpha = pow(alpha_liquid_vapor(t, HDO, CondensedOverVapor), sp(2.0));
        break;
      case H216O:
      default:
        // alpha stays 1 (non-fractionating).
        break;
    }

    // KINETIC HOOK: a future kinetic evaporation factor (iCAM wiso_akel) would
    // modify `alpha` here before the direction is applied.
    return (dir == VaporOverCondensed) ? (sp(1) / alpha) : alpha;
  }

  // -----------------------------------------------------------------------
  // Ice(solid) <-> vapor equilibrium fractionation factor.
  //
  // Merlivat & Nief (1967) for HDO; Majoube (1971b) for H2-18O. Single
  // functional form for all species; T in Kelvin:
  //   alpha = exp( a/T^2 + b/T + c )
  // Coefficients transcribed directly from water_isotopes.F90 (lines 183-185).
  //
  // Note: in iCAM this equilibrium factor is applied on DEPOSITION
  // (vapor->ice); SUBLIMATION (ice->vapor) is treated as non-fractionating
  // (alpha = 1) by the caller, not here. This function is the pure equilibrium
  // factor; the deposition/sublimation choice is process logic.
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_ice_vapor(const ScalarT& t,
                                 const WisoSpecies species,
                                 const WisoAlphaDir dir)
  {
    // alpha = exp(a/T^2 + b/T + c)
    //                       a               b              c
    constexpr Real hdo_a = sp(16289.0);
    constexpr Real hdo_b = sp(0.0);
    constexpr Real hdo_c = sp(-9.45e-2);

    constexpr Real o18_a = sp(0.0);
    constexpr Real o18_b = sp(11.839);
    constexpr Real o18_c = sp(-28.224e-3);

    ScalarT alpha(1);

    switch (species) {
      case HDO: {
        const ScalarT it  = sp(1) / t;
        const ScalarT it2 = it * it;
        alpha = exp(hdo_a * it2 + hdo_b * it + hdo_c);
        break;
      }
      case H218O: {
        const ScalarT it  = sp(1) / t;
        const ScalarT it2 = it * it;
        alpha = exp(o18_a * it2 + o18_b * it + o18_c);
        break;
      }
      case H217O:
        // Mass-dependent from H2-18O (Schoenemann et al. 2014).
        alpha = pow(alpha_ice_vapor(t, H218O, CondensedOverVapor), sp(0.529));
        break;
      case HTO:
        // From HDO (isoCAM3).
        alpha = pow(alpha_ice_vapor(t, HDO, CondensedOverVapor), sp(2.0));
        break;
      case H216O:
      default:
        // alpha stays 1 (non-fractionating).
        break;
    }

    // KINETIC HOOK: a future ice-condensation factor (iCAM wiso_akci, with the
    // supersaturation function and T < ~253 K mask) would modify `alpha` here
    // before the direction is applied.
    return (dir == VaporOverCondensed) ? (sp(1) / alpha) : alpha;
  }
};

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP
