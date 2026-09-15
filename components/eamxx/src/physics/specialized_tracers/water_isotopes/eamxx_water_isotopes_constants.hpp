#ifndef EAMXX_WATER_ISOTOPES_CONSTANTS_HPP
#define EAMXX_WATER_ISOTOPES_CONSTANTS_HPP

#include "share/core/eamxx_types.hpp"  // for scream::Real and scream::sp()

namespace scream {
namespace wiso {

/*
 * Physical constants for water isotopologues.
 *
 * Ported from the Fortran module water_isotopes.F90
 * These constants define isotopic fractionation behavior, molecular properties,
 * and reference ratios for water isotope tracers.
 *
 * All species-specific constants are arrays indexed by WisoSpecies enum:
 *   H216O = 0 ("ordinary water," non-fractionating)
 *   HDO   = 1 (singly deuterated water, HD16O)
 *   H218O = 2 (oxygen-18 substituted water, H218O)
 *   H217O = 3 (oxygen-17 substituted water, H217O)
 *   HTO   = 4 (tritiated water, HT16O)
 *
 * Formulations are now selected at runtime via WaterIsotopeRuntimeOptions.
 * See README.md for details on available formulations and their scientific references.
 *
 * Default formulations: Horita & Wesolowski 1994 (liquid/vapor), Merlivat 1978 (diffusivity),
 *                       Normalized ratios, No ocean enrichment,
 *                       Merlivat & Nief 1967 + Majoube 1971b (ice/vapor)
 */

// ============================================================================
// Runtime configuration enums and struct
// ============================================================================

// Liquid/vapor equilibrium fractionation formulation
enum class LiquidVaporFractionation {
  HoritaWesolowski1994 = 0,  // Default: Horita & Wesolowski (1994)
  Majoube1971 = 1,           // Alternative: Majoube (1971)
  FormulationCount
};

// Diffusivity ratio formulation
enum class DiffusivityFormulation {
  Merlivat1978 = 0,  // Default: Merlivat (1978)
  Cappa2003 = 1,     // Alternative: Cappa et al. (2003)
  FormulationCount
};

// Standard isotope ratio formulation
enum class StandardRatioFormulation {
  Normalized = 0,        // Default: All 1.0 for best numerics
  NaturalAbundance = 1,  // Alternative: Natural abundance values
  FormulationCount
};

// Ocean surface enrichment formulation
enum class OceanEnrichmentFormulation {
  None = 0,  // Default: No enrichment (all 1.0)
  LGM = 1,   // Alternative: Last Glacial Maximum values
  FormulationCount
};

// Ice/vapor equilibrium fractionation formulation
enum class IceVaporFractionation {
  MerlivatNief1967 = 0,  // Default: Merlivat & Nief (1967) HDO + Majoube (1971) O18
  IsoCAM3 = 1,           // Alternative: isoCAM3 formulation
  FormulationCount
};

// Runtime configuration struct - holds user's formulation choices
struct WaterIsotopeRuntimeOptions {
  LiquidVaporFractionation liquid_vapor = LiquidVaporFractionation::HoritaWesolowski1994;
  DiffusivityFormulation diffusivity = DiffusivityFormulation::Merlivat1978;
  StandardRatioFormulation standard_ratio = StandardRatioFormulation::Normalized;
  OceanEnrichmentFormulation ocean_enrichment = OceanEnrichmentFormulation::None;
  IceVaporFractionation ice_vapor = IceVaporFractionation::MerlivatNief1967;
};

template <typename Scalar>
struct WaterIsotopeConstants
{
  using Real = Scalar;

  // Number of isotope species
  static constexpr int num_species = 5;

  // -----------------------------------------------------------------------
  // Active constants (runtime-selected, non-constexpr)
  // -----------------------------------------------------------------------

  // Diffusivity ratios (D_isotope/D_H2O)
  Scalar difrm[num_species];

  // Model standard isotope ratios
  Scalar rstd[num_species];

  // Ocean surface enrichment
  Scalar boce[num_species];

  // Liquid/vapor equilibrium fractionation coefficients - requires 5 coefficients
  // Coefficients for alpha = exp(polynomial in T)
  // HDO uses: alpha = exp(a*T^3 + b*T^2 + c*T + d + e/T^3)
  // H218O uses: alpha = exp(a/T^3 + b/T^2 + c/T + d)
  // H217O and HTO are derived from these by power laws
  // Not all formulations have all coefficients, depends on regression.
  static constexpr Scalar AlphaLiqVap_CoefA_table
      [LiquidVaporFractionation::FormulationCount][num_species] = {
        { 0.0, 1158.8e-12, 0.35041e6, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,   24.844e3,   1.137e3, 0.0, 0.0 }  // Majoube 1971
  };
  static constexpr Scalar AlphaLiqVap_CoefB_table
      [LiquidVaporFractionation::FormulationCount][num_species] = {
        { 0.0, -1620.1e-9, -1.6664e3, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,    -76.248,   -0.4156, 0.0, 0.0 }  // Majoube 1971
  };
  static constexpr Scalar AlphaLiqVap_CoefC_table
      [LiquidVaporFractionation::FormulationCount][num_species] = {
        { 0.0,  794.84e-6,    6.7123, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,  52.612e-3,-2.0667e-3, 0.0, 0.0 }  // Majoube 1971
  };
  static constexpr Scalar AlphaLiqVap_CoefD_table
      [LiquidVaporFractionation::FormulationCount][num_species] = {
        { 0.0, -161.04e-3, -7.685e-3, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,        0.0,       0.0, 0.0, 0.0 }  // Majoube 1971
  };
  static constexpr Scalar AlphaLiqVap_CoefE_table
      [LiquidVaporFractionation::FormulationCount][num_species] = {
        { 0.0,   2.9992e6,       0.0, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,        0.0,       0.0, 0.0, 0.0 }  // Majoube 1971
  };

  // Ice/vapor equilibrium fractionation coefficients
  // Coefficients for alpha = exp(a/T^2 + b/T + c)
  Scalar alpai[num_species];  // Coefficient A
  Scalar alpbi[num_species];  // Coefficient B
  Scalar alpci[num_species];  // Coefficient C

  static constexpr Scalar AlphaIceVap_CoefA_table
      [IceVaporFractionation::FormulationCount][num_species] = {
        { 0.0, 16289.0, 0.0, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,   16288.0,   0.0, 0.0, 0.0 }  // Majoube 1971
  };  
  static constexpr Scalar AlphaIceVap_CoefB_table
      [IceVaporFractionation::FormulationCount][num_species] = {
        { 0.0, 0.0, 11.839, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,   0.0,   11.839, 0.0, 0.0 }  // Majoube 1971
  };  
  static constexpr Scalar AlphaIceVap_CoefC_table
      [IceVaporFractionation::FormulationCount][num_species] = {
        { 0.0, -9.45e-2, -28.224e-3, 0.0, 0.0 }, // Horita and Wesolowski 1994
        { 0.0,  -9.34e-2,   -28.224e-3, 0.0, 0.0 }  // Majoube 1971
  };


  // -----------------------------------------------------------------------
  // Constructors
  // -----------------------------------------------------------------------

  // Constructor with runtime options (selects formulations)
  KOKKOS_INLINE_FUNCTION
  WaterIsotopeConstants(const WaterIsotopeRuntimeOptions& opts) {
    // Select diffusivity formulation
    const Scalar* diff_src = (opts.diffusivity == DiffusivityFormulation::Cappa2003)
                              ? difrm_cappa2003 : difrm_merlivat1978;
    copy_array(difrm, diff_src, num_species);

    // Select standard ratio formulation
    const Scalar* ratio_src = (opts.standard_ratio == StandardRatioFormulation::NaturalAbundance)
                               ? rstd_natural_abundance : rstd_normalized;
    copy_array(rstd, ratio_src, num_species);

    // Select ocean enrichment formulation
    const Scalar* ocean_src = (opts.ocean_enrichment == OceanEnrichmentFormulation::LGM)
                               ? boce_lgm : boce_none;
    copy_array(boce, ocean_src, num_species);

    // Select coefficients for liquid/vapor formulation selected
    // TODO: rename apla*l variables throughout.
    KOKKOS_INLINE_FUNCTION Scalar alpal(int s) const { return AlphaLiqVap_CoefA_table[int(opts.liquid_vapor)][s] } 
    KOKKOS_INLINE_FUNCTION Scalar alpbl(int s) const { return AlphaLiqVap_CoefB_table[int(opts.liquid_vapor)][s] }
    KOKKOS_INLINE_FUNCTION Scalar alpcl(int s) const { return AlphaLiqVap_CoefC_table[int(opts.liquid_vapor)][s] }
    KOKKOS_INLINE_FUNCTION Scalar alpdl(int s) const { return AlphaLiqVap_CoefD_table[int(opts.liquid_vapor)][s] }
    KOKKOS_INLINE_FUNCTION Scalar alpel(int s) const { return AlphaLiqVap_CoefE_table[int(opts.liquid_vapor)][s] }
    // Ice/vapor
    KOKKOS_INLINE_FUNCTION Scalar alpai(int s) const { return AlphaIceVap_CoefA_table[int(opts.ice_vapor)][s] } 
    KOKKOS_INLINE_FUNCTION Scalar alpbi(int s) const { return AlphaIceVap_CoefB_table[int(opts.ice_vapor)][s] }
    KOKKOS_INLINE_FUNCTION Scalar alpci(int s) const { return AlphaIceVap_CoefC_table[int(opts.ice_vapor)][s] }
  }

  // Default constructor (uses default formulations)
  KOKKOS_INLINE_FUNCTION
  WaterIsotopeConstants() : WaterIsotopeConstants(WaterIsotopeRuntimeOptions{}) {}

  // -----------------------------------------------------------------------
  // Molecular properties (species-indexed arrays) - always same
  // -----------------------------------------------------------------------

  // Isotopic substitutions (mass-dependent factor)
  // Ported from CAM6 water_isotopes.F90
  // C++: expanded to 5 elements with H217O and HTO derived values
  static constexpr Scalar fisub[num_species] = {
    1.0,  // H216O (non-fractionating)
    2.0,  // HDO (deuterium substitution)
    1.0,  // H218O (oxygen substitution)
    1.0,  // H217O (oxygen substitution, same as H218O)
    3.0   // HTO (tritium substitution)
  };

private:
  // -----------------------------------------------------------------------
  // All formulation data as static constexpr (compile-time constants)
  // -----------------------------------------------------------------------

  // Diffusivity ratios - Merlivat 1978 (default)
  static constexpr Scalar difrm_merlivat1978[num_species] = {
    1.0,     // H216O
    0.9757,  // HDO
    0.9727,  // H218O
    0.9727,  // H217O (TODO: assumed same as H218O)
    0.9757   // HTO (TODO: assumed same as HDO)
  };

  // Diffusivity ratios - Cappa et al. 2003
  static constexpr Scalar difrm_cappa2003[num_species] = {
    1.0,     // H216O
    0.9839,  // HDO
    0.9691,  // H218O
    0.9691,  // H217O (TODO: assumed same as H218O)
    0.9839   // HTO (TODO: assumed same as HDO)
  };

  // Standard ratios - Normalized (default)
  static constexpr Scalar rstd_normalized[num_species] = {
    1.0,  // H216O
    1.0,  // HDO
    1.0,  // H218O
    1.0,  // H217O
    1.0   // HTO
  };

  // Standard ratios - Natural abundance
  static constexpr Scalar rstd_natural_abundance[num_species] = {
    1.0,         // H216O (reference)
    155.76e-6,   // HDO
    2005.20e-6,  // H218O
    379.9e-6,    // H217O
    1.0          // HTO (TODO: define more robustly)
  };

  // Ocean surface enrichment - None (default)
  static constexpr Scalar boce_none[num_species] = {
    1.0,  // H216O
    1.0,  // HDO
    1.0,  // H218O
    1.0,  // H217O
    1.0   // HTO
  };

  // Ocean surface enrichment - LGM
  static constexpr Scalar boce_lgm[num_species] = {
    1.0,      // H216O
    1.0128,   // HDO
    1.0016,   // H218O
    1.0008,   // H217O (TODO)
    1.0       // HTO (TODO)
  };

  // -----------------------------------------------------------------------
  // Ice/vapor equilibrium fractionation coefficients
  // -----------------------------------------------------------------------

  
  // Helper function to copy arrays
  KOKKOS_INLINE_FUNCTION
  void copy_array(Scalar* dest, const Scalar* src, int n) {
    for (int i = 0; i < n; ++i) {
      dest[i] = src[i];
    }
  }

public:

  // -----------------------------------------------------------------------
  // Kinetic fractionation parameters (Merlivat & Jouzel method)
  // These do not vary by formulation
  // -----------------------------------------------------------------------

  // Surface kinetic exchange
  // From water_isotopes.F90
  // Note: H216O entries are 0.0 (non-fractionating)
  static constexpr Scalar aksmc[num_species] = {
    0.0,      // H216O
    0.00528,  // HDO
    0.006,    // H218O
    0.006,    // H217O (assumed same as H218O)
    0.00528   // HTO (assumed same as HDO)
  };

  // Parameter A for kinetic fractionation
  static constexpr Scalar akrfa[num_species] = {
    0.0,        // H216O
    0.2508e-3,  // HDO
    0.285e-3,   // H218O
    0.285e-3,   // H217O (assumed same as H218O)
    0.2508e-3   // HTO (assumed same as HDO)
  };

  // Parameter B for kinetic fractionation
  static constexpr Scalar akrfb[num_species] = {
    0.0,        // H216O
    0.7216e-3,  // HDO
    0.82e-3,    // H218O
    0.82e-3,    // H217O (assumed same as H218O)
    0.7216e-3   // HTO (assumed same as HDO)
  };

  // -----------------------------------------------------------------------
  // Physical constants for kinetic calculations
  // -----------------------------------------------------------------------
  // RPF - check if these already exist in physics constants, and if 
  // there are more precise values available.
  // Molecular diffusivity of air [m2/s]
  // From water_isotopes.F90 line 215
  static constexpr Scalar MolecularDiffusivityAir = 2.36e-5;

  // Dynamic viscosity of air [Pa*s or kg/(m*s)]
  // From water_isotopes.F90 line 216
  // Assumed constant, TODO replace with temperature-dependent formulation
  static constexpr Scalar DynamicViscosityAir = 1.7e-5;
};

// Convenience alias for Real precision
using WaterIsotopeConstantsReal = WaterIsotopeConstants<Real>;

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_CONSTANTS_HPP
