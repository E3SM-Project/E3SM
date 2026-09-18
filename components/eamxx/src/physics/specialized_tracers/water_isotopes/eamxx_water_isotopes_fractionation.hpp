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
 * Derived from the iCAM Fortran module water_isotopes.F90 (functions wiso_alpl
 * and wiso_alpi, original author David Noone). These are pure, device-callable functions
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
 */

// Water isotopologues. HDO and H218O are computed directly from the tables;
// H217O and HTO are derived by mass-dependent power laws. H216O (ordinary
// water) is non-fractionating (alpha == 1).
enum WisoSpecies {
  H216O = 0,  // ordinary water; alpha == 1
  HDO   = 1,  // HD16O (deuterium)
  H218O = 2,  // H218O (oxygen-18)
  H217O = 3,  // H217O; = alpha(H218O)^0.529 (Schoenemann et al. 2014)
  HTO   = 4   // HT16O (tritiated water); = alpha(HDO)^2.0 (isoCAM3 assumption)
};

// Which R-ratio the returned factor represents.
enum WisoAlphaDir {
  CondensedOverVapor = 0,  // raw table value R_condensed/R_vapor (>= 1)
  VaporOverCondensed = 1   // reciprocal, R_vapor/R_condensed (<= 1)
};

struct WaterIsotopeFractionation
{
// Private helper function to calculate alpha_eq given species, temperature
private:
  // Mass-dependent scaling exponents
  static constexpr double H217O_exponent = 0.529;  // Schoenemann et al. (2014)
  static constexpr double HTO_exponent = 2.0;      // isoCAM3 assumption

  // Common fractionation logic: derived species, direction handling
  template <typename ScalarT, typename BaseFunc>
  KOKKOS_INLINE_FUNCTION
  static ScalarT compute_alpha(const ScalarT& t,
                                const WisoSpecies species,
                                const WisoAlphaDir dir,
                                BaseFunc base_alpha)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    ScalarT alpha(1);

    switch (species) {
      case HDO:
        alpha = base_alpha(t, HDO);
        break;
      case H218O:
        alpha = base_alpha(t, H218O);
        break;
      case H217O:
        // Derived from H218O via mass-dependent fractionation
        alpha = pow(base_alpha(t, H218O), RealT(H217O_exponent));
        break;
      case HTO:
        // Derived from HDO via mass-dependent fractionation
        alpha = pow(base_alpha(t, HDO), RealT(HTO_exponent));
        break;
      case H216O:
      default:
        // Non-fractionating (alpha = 1)
        break;
    }

    // Apply direction
    return (dir == VaporOverCondensed) ? (RealT(1) / alpha) : alpha;
  }
public:
  // -----------------------------------------------------------------------
  // Liquid <-> vapor equilibrium fractionation factor (with runtime constants).
  //
  // Two functional forms by species; T in Kelvin.
  // Coefficients come from the constants struct (selected formulation).
  // -----------------------------------------------------------------------
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_liquid_vapor(
      const ScalarT& t,
      const WisoSpecies species,
      const WisoAlphaDir dir,
      const WaterIsotopeConstants<typename ekat::ScalarTraits<ScalarT>::scalar_type>& constants)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    // Define the liquid-vapor polynomial (the unique part)
    auto base = [&](const ScalarT& temp, WisoSpecies sp) -> ScalarT {
      if (sp == HDO) {
        // HDO: alpha = exp(a*T³ + b*T² + c*T + d + e/T³)
        const ScalarT t2 = temp * temp;
        const ScalarT t3 = t2 * temp;
        return exp(constants.alpal(HDO) * t3 +
                   constants.alpbl(HDO) * t2 +
                   constants.alpcl(HDO) * temp +
                   constants.alpdl(HDO) +
                   constants.alpel(HDO) / t3);
      } else {  // H218O
        // H218O: alpha = exp(a/T³ + b/T² + c/T + d)
        const ScalarT it  = RealT(1) / temp;
        const ScalarT it2 = it * it;
        const ScalarT it3 = it2 * it;
        return exp(constants.alpal(H218O) * it3 +
                   constants.alpbl(H218O) * it2 +
                   constants.alpcl(H218O) * it +
                   constants.alpdl(H218O));
      }
    };
  
    // Use the common helper
    return compute_alpha(t, species, dir, base);
  }

  // Backward-compatible overload using default formulation (Horita & Wesolowski 1994)
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_liquid_vapor(const ScalarT& t,
                                    const WisoSpecies species,
                                    const WisoAlphaDir dir)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    WaterIsotopeConstants<RealT> constants;  // Uses default formulations
    return alpha_liquid_vapor(t, species, dir, constants);
  }

  // -----------------------------------------------------------------------
  // Ice(solid) <-> vapor equilibrium fractionation factor (with runtime constants).
  //
  // Single functional form for all species; T in Kelvin:
  //   alpha = exp( a/T^2 + b/T + c )
  // Coefficients come from the constants struct (selected formulation).
  //
  // -----------------------------------------------------------------------
    template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_ice_vapor(
      const ScalarT& t,
      const WisoSpecies species,
      const WisoAlphaDir dir,
      const WaterIsotopeConstants<typename ekat::ScalarTraits<ScalarT>::scalar_type>& constants)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;

    // Define the ice-vapor polynomial (the unique part)
    auto base = [&](const ScalarT& temp, WisoSpecies sp) -> ScalarT {
      // Both species use same form: alpha = exp(a/T² + b/T + c)
      const ScalarT it  = RealT(1) / temp;
      const ScalarT it2 = it * it;
      return exp(constants.alpai(sp) * it2 +
                 constants.alpbi(sp) * it +
                 constants.alpci(sp));
    };

    // Use the common helper
    return compute_alpha(t, species, dir, base);
  }

  // Backward-compatible overload using default formulation
  template <typename ScalarT>
  KOKKOS_INLINE_FUNCTION
  static ScalarT alpha_ice_vapor(const ScalarT& t,
                                 const WisoSpecies species,
                                 const WisoAlphaDir dir)
  {
    using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
    WaterIsotopeConstants<RealT> constants;  // Uses default formulations
    return alpha_ice_vapor(t, species, dir, constants);
  }
};

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_FRACTIONATION_HPP
