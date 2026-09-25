#ifndef EAMXX_WATER_ISOTOPES_CONSTANTS_HPP
#define EAMXX_WATER_ISOTOPES_CONSTANTS_HPP

#include "share/core/eamxx_types.hpp"  // for scream::Real
#include "share/util/eamxx_utils.hpp"

namespace scream {
namespace wiso {

/*
 * Physical constants for water isotopologues.
 *
 * Ported from the Fortran module water_isotopes.F90
 * These constants define isotopic fractionation behavior, molecular properties,
 * and reference ratios for water isotope tracers.
 *
 * All species-specific constants are arrays indexed by WaterIsotopologues enum:
 *   H216O = 0 ("ordinary water," non-fractionating)
 *   HDO   = 1 (singly deuterated water, HD16O)
 *   H218O = 2 (oxygen-18 substituted water, H218O)
 *   H217O = 3 (oxygen-17 substituted water, H217O)
 *   HTO   = 4 (tritiated water, HT16O)
 *
 * Formulations are now selected at runtime via WaterIsotopeRuntimeOptions.
 * See README.md for details on available formulations and their scientific references.
 *
 * Default formulations: Horita & Wesolowski 1994 (liquid/vapor),
 *                       Normalized ratios, No ocean enrichment,
 *                       Merlivat & Nief 1967 + Majoube 1971 (ice/vapor)
 */

// Define water isotopologue species
enum class WaterIsotopologues {
  H216O = 0,  // ordinary water; alpha == 1
  HDO   = 1,  // HD16O (deuterium)
  H218O = 2,  // H218O (oxygen-18)
  H217O = 3,  // H217O; = alpha(H218O)^0.529 (Schoenemann et al. 2014)
  HTO   = 4,  // HT16O (tritiated water); = alpha(HDO)^2.0 (isoCAM3 assumption)
  Count
};

// Substituted element. Indexes the coefficient tables below.
enum class IsoElement { Hydrogen = 0, Oxygen = 1, Count };

// Condensed phase of the vapor <-> condensate equilibrium.
enum class CondensedPhase { Liquid = 0, Ice = 1, Count };

/* Coefficients of the temperature polynomial:

     10^3 * ln(alpha) = T3*T^3 + T2*T^2 + T1*T + T0
                        + T_1/T + T_2/T^2 + T_3/T^3 + T_4/T^4 + T_6/T^6

   Not every formulation uses every term; unused terms are 0. T_4 and T_6 are
   currently zero in every formulation below, but are used in coefficient
   sets to be added. */
struct PolynomialCoefficients {
  Real T3, T2, T1, T0, T_1, T_2, T_3, T_4, T_6;
};

// set boundaries from polynomial regression
struct TemperatureBounds {
  Real Tmin, Tmax;  // [K] range of the published regression
};

struct EquilibriumFractionationCoefficients {
  TemperatureBounds      tbounds;
  PolynomialCoefficients coeffs;
};

// Liquid/vapor equilibrium fractionation formulation options
enum class LiquidVaporFractionation {
  HoritaWesolowski1994 = 0,  // Default: Horita & Wesolowski (1994)
  Majoube1971 = 1,           // Alternative: Majoube (1971)
  FormulationCount
};

// Ice/vapor equilibrium fractionation formulation options
enum class IceVaporFractionation {
  MerlivatNief1967 = 0,  // Default: Merlivat & Nief (1967) HDO + Majoube (1971) O18
  IsoCAM3 = 1,           // Alternative: isoCAM3 formulation
  FormulationCount
};

// Coefficient tables, indexed [formulation][element]. 
static constexpr EquilibriumFractionationCoefficients
alpha_eq_liq_table[etoi(LiquidVaporFractionation::FormulationCount)]
                  [etoi(IsoElement::Count)] = {
  // Horita & Wesolowski (1994), 0-364 C
  { /* Hydrogen */ {{273.15, 637.15}, {1.1588e-6, -1.6201e-3, 7.9484e-1, -1.6104e2, 0., 0., 2.9992e9, 0., 0.}},
    /* Oxygen   */ {{273.15, 637.15}, {0., 0., 0., -7.685, 6.7123e3, -1.6664e6, 3.5041e8, 0., 0.}} },
  // Majoube (1971), 0-100 C
  { /* Hydrogen */ {{273.15, 373.15}, {0., 0., 0., 5.2612e1, -7.6248e4, 2.4844e7, 0., 0., 0.}},
    /* Oxygen   */ {{273.15, 373.15}, {0., 0., 0., -2.0667, -4.156e2, 1.137e6, 0., 0., 0.}} }
};

static constexpr EquilibriumFractionationCoefficients
alpha_eq_ice_table[etoi(IceVaporFractionation::FormulationCount)]
                  [etoi(IsoElement::Count)] = {
  // Default: two different studies, one per element.
  { /* Hydrogen */ {{233.15, 273.15}, {0., 0., 0., -9.45e1, 0., 1.6289e7, 0., 0., 0.}},  // Merlivat & Nief (1967), -40..0 C
    /* Oxygen   */ {{239.75, 273.15}, {0., 0., 0., -2.8224e1, 1.1839e4, 0., 0., 0., 0.}} },  // Majoube (1971), NOT M&N
  // isoCAM3: same functional fits, extrapolated to -70 C.
  { /* Hydrogen */ {{203.15, 273.15}, {0., 0., 0., -9.34e1, 0., 1.6288e7, 0., 0., 0.}},
    /* Oxygen   */ {{203.15, 273.15}, {0., 0., 0., -2.8224e1, 1.1839e4, 0., 0., 0., 0.}} }
};

// Provenance strings, for logging. Host-only: deliberately not part of the
// device-side payload (carrying them would grow the kernel closure).
static constexpr const char* liquid_vapor_ref
    [etoi(LiquidVaporFractionation::FormulationCount)] = {
  "Horita & Wesolowski (1994)",
  "Majoube (1971)"
};
static constexpr const char* ice_vapor_ref
    [etoi(IceVaporFractionation::FormulationCount)] = {
  "Merlivat & Nief (1967) [HDO] + Majoube (1971) [H218O]",
  "isoCAM3"
};

// ============================================================================
// Runtime configuration enums and struct
// ============================================================================

// Standard isotope ratio formulation
enum class StandardRatioFormulation {
  Normalized = 0,        // Default: All 1.0 for best numerics
  NaturalAbundance = 1,  // Alternative: Natural abundance values
  FormulationCount
};

// Ocean surface enrichment formulation
enum class OceanEnrichmentFormulation {
  Modern = 0,// Default: No enrichment (all 1.0)
  LGM = 1,   // Alternative: Last Glacial Maximum values
  FormulationCount
};

// Runtime configuration struct - holds user's formulation choices
struct WaterIsotopeRuntimeOptions {
  LiquidVaporFractionation liquid_vapor = LiquidVaporFractionation::HoritaWesolowski1994;
  StandardRatioFormulation standard_ratio = StandardRatioFormulation::Normalized;
  OceanEnrichmentFormulation ocean_enrichment = OceanEnrichmentFormulation::Modern;
  IceVaporFractionation ice_vapor = IceVaporFractionation::MerlivatNief1967;
};

template <typename Scalar>
struct WaterIsotopeConstants
{
  using Real = Scalar;

  // Number of isotope species
  static constexpr int num_species = etoi(WaterIsotopologues::Count);

  // -----------------------------------------------------------------------
  // Active constants (runtime-selected)
  // -----------------------------------------------------------------------

  // Model standard isotope ratios
  static constexpr Scalar rstd_table
      [etoi(StandardRatioFormulation::FormulationCount)][num_species] = {
        { 1.0,       1.0,        1.0,      1.0, 1.0 }, // Normalized by VSMOW
        { 1.0, 155.76e-6, 2005.20e-6, 379.9e-6, 1.0 } // natural abundance
  };

  // Ocean surface enrichment
  static constexpr Scalar boce_table
    [etoi(OceanEnrichmentFormulation::FormulationCount)][num_species] = {
      { 1.0,    1.0,    1.0,    1.0, 1.0 }, // Modern VSMOW
      { 1.0, 1.0128, 1.0016, 1.0008, 1.0} // LGM
  };

  // -----------------------------------------------------------------------
  // Runtime-selected formulation options
  // -----------------------------------------------------------------------

  WaterIsotopeRuntimeOptions opts_;

private:
  // Equilibrium fractionation coefficients for the *selected* formulations,
  // resolved once on construction and stored by value.
  EquilibriumFractionationCoefficients
    eq_[etoi(CondensedPhase::Count)][etoi(IsoElement::Count)];

public:
  // -----------------------------------------------------------------------
  // Constructors
  // -----------------------------------------------------------------------

  // Constructor with runtime options (selects formulations)
  KOKKOS_INLINE_FUNCTION
  WaterIsotopeConstants(const WaterIsotopeRuntimeOptions& opts) : opts_(opts) {
    for (int e = 0; e < etoi(IsoElement::Count); ++e) {
      eq_[etoi(CondensedPhase::Liquid)][e] =
          alpha_eq_liq_table[etoi(opts.liquid_vapor)][e];
      eq_[etoi(CondensedPhase::Ice)][e] =
          alpha_eq_ice_table[etoi(opts.ice_vapor)][e];
    }
  }

  // Default constructor (uses default formulations)
  KOKKOS_INLINE_FUNCTION
  WaterIsotopeConstants() : WaterIsotopeConstants(WaterIsotopeRuntimeOptions{}) {}

  // -----------------------------------------------------------------------
  // Accessors (select table row based on stored runtime options)
  // -----------------------------------------------------------------------

  // Select standard ratio formulation
  KOKKOS_INLINE_FUNCTION
  Scalar ratio_src(WaterIsotopologues s) const { return rstd_table[int(opts_.standard_ratio)][etoi(s)]; }

  // Select ocean enrichment formulation
  KOKKOS_INLINE_FUNCTION
  Scalar ocean_src(WaterIsotopologues s) const { return boce_table[int(opts_.ocean_enrichment)][etoi(s)]; }

  // Which element's coefficient set a species uses.
  KOKKOS_INLINE_FUNCTION
  static IsoElement element_of(WaterIsotopologues s) {
    return (s == WaterIsotopologues::HDO || s == WaterIsotopologues::HTO) ? IsoElement::Hydrogen : IsoElement::Oxygen;
  }

  // Equilibrium fractionation polynomial coefficients 
  KOKKOS_INLINE_FUNCTION
  const PolynomialCoefficients& alpha_eq_coeffs(CondensedPhase p, IsoElement e) const {
    return eq_[etoi(p)][etoi(e)].coeffs;
  }

  // Temperature range [K] over which the selected regression was fitted.
  // Outside it the polynomial is an extrapolation, not a fit.
  KOKKOS_INLINE_FUNCTION
  const TemperatureBounds& tbounds(CondensedPhase p, IsoElement e) const {
    return eq_[etoi(p)][etoi(e)].tbounds;
  }

};

// Convenience alias for Real precision
using WaterIsotopeConstantsReal = WaterIsotopeConstants<Real>;

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_CONSTANTS_HPP
