#ifndef EAMXX_WATER_ISOTOPES_CONSTANTS_HPP
#define EAMXX_WATER_ISOTOPES_CONSTANTS_HPP

#include "share/core/eamxx_types.hpp"  // for scream::Real and scream::sp()

namespace scream {
namespace wiso {

/*
 * Physical constants for water isotopologues.
 *
 * Ported from the Fortran module water_isotopes.F90 (share/util/water_isotopes.F90).
 * These constants define isotopic fractionation behavior, molecular properties,
 * and reference ratios for water isotope tracers.
 *
 * All species-specific constants are arrays indexed by WisoSpecies enum:
 *   H216O = 0 (ordinary water, non-fractionating)
 *   HDO   = 1 (deuterated water, HD16O)
 *   H218O = 2 (oxygen-18 water, H2-18O)
 *   H217O = 3 (oxygen-17 water, H2-17O) - derived from H218O
 *   HTO   = 4 (tritiated water, HT16O) - derived from HDO
 *
 * Alternative formulations can be selected via preprocessor macros:
 *   WISO_USE_MAJOUBE_1971     - Use Majoube 1971a liquid/vapor coefficients
 *   WISO_USE_CAPPA_2003       - Use Cappa et al. 2003 diffusivity ratios
 *   WISO_USE_NATURAL_ABUNDANCE - Use natural abundance standard ratios
 *   WISO_USE_LGM_OCEAN        - Use LGM ocean enrichment values
 *   WISO_USE_ISOCAM3_ICE      - Use isoCAM3 ice/vapor coefficients
 *
 * Default (no macros): Horita & Wesolowski 1994, Merlivat 1978 diffusivity,
 *                      Merlivat & Nief 1967 / Majoube 1971b ice/vapor
 */

template <typename Scalar>
struct WaterIsotopeConstants
{
  using Real = Scalar;

  // Number of isotope species
  static constexpr int num_species = 5;

  // -----------------------------------------------------------------------
  // Molecular properties (species-indexed arrays)
  // -----------------------------------------------------------------------

  // Isotopic substitutions (mass-dependent factor)
  // From water_isotopes.F90 lines 114-115
  // Fortran: fisub = (/ 1.0, 1.0, 2.0, 1.0 /) for H2O, H216O, HDO, H218O
  // C++: expanded to 5 elements with H217O and HTO derived values
  static constexpr Scalar fisub[num_species] = {
    sp(1.0),  // H216O (non-fractionating)
    sp(2.0),  // HDO (deuterium substitution)
    sp(1.0),  // H218O (oxygen substitution)
    sp(1.0),  // H217O (oxygen substitution, same as H218O)
    sp(3.0)   // HTO (tritium substitution)
  };

  // Diffusivity ratios (D_isotope/D_H2O)
  // From water_isotopes.F90 lines 119-125
  // Default: Merlivat 1978 direct from paper (line 124)
  // Alternative: WISO_USE_CAPPA_2003 for Cappa et al. 2003 values
#ifdef WISO_USE_CAPPA_2003
  // Cappa et al. 2003 (Fortran line 125)
  static constexpr Scalar difrm[num_species] = {
    sp(1.0),     // H216O
    sp(0.9839),  // HDO
    sp(0.9691),  // H218O
    sp(0.9691),  // H217O (assumed same as H218O)
    sp(0.9839)   // HTO (assumed same as HDO)
  };
#else
  // Merlivat 1978 direct from paper (default, Fortran line 124)
  static constexpr Scalar difrm[num_species] = {
    sp(1.0),     // H216O
    sp(0.9757),  // HDO
    sp(0.9727),  // H218O
    sp(0.9727),  // H217O (assumed same as H218O)
    sp(0.9757)   // HTO (assumed same as HDO)
  };
#endif

  // -----------------------------------------------------------------------
  // Standard ratios (species-indexed arrays)
  // -----------------------------------------------------------------------

  // Model standard isotope ratio
  // From water_isotopes.F90 lines 128-134
  // Default: all 1.0 for best numerics (Fortran line 130)
  // Alternative: WISO_USE_NATURAL_ABUNDANCE for natural abundance values
#ifdef WISO_USE_NATURAL_ABUNDANCE
  // Natural abundance (Fortran line 132)
  static constexpr Scalar rstd[num_species] = {
    sp(1.0),         // H216O (reference)
    sp(155.76e-6),   // HDO
    sp(2005.20e-6),  // H218O
    sp(379.9e-6),    // H217O (approximately 0.189 * H218O)
    sp(1.0e-18)      // HTO (trace amount)
  };
#else
  // All 1.0 for best numerics (default, Fortran line 130)
  static constexpr Scalar rstd[num_species] = {
    sp(1.0),  // H216O
    sp(1.0),  // HDO
    sp(1.0),  // H218O
    sp(1.0),  // H217O
    sp(1.0)   // HTO
  };
#endif

  // Ocean surface enrichment
  // From water_isotopes.F90 lines 136-140
  // Default: all 1.0 (no enrichment, Fortran line 140)
  // Alternative: WISO_USE_LGM_OCEAN for Last Glacial Maximum values
#ifdef WISO_USE_LGM_OCEAN
  // LGM values (Fortran line 139)
  static constexpr Scalar boce[num_species] = {
    sp(1.0),      // H216O
    sp(1.0128),   // HDO
    sp(1.0016),   // H218O
    sp(1.0008),   // H217O (interpolated)
    sp(1.00671)   // HTO (assumed)
  };
#else
  // No enrichment (default, Fortran line 140)
  static constexpr Scalar boce[num_species] = {
    sp(1.0),  // H216O
    sp(1.0),  // HDO
    sp(1.0),  // H218O
    sp(1.0),  // H217O
    sp(1.0)   // HTO
  };
#endif

  // -----------------------------------------------------------------------
  // Kinetic fractionation parameters (Merlivat & Jouzel method)
  // -----------------------------------------------------------------------

  // Surface kinetic exchange
  // From water_isotopes.F90 lines 144-147
  // Note: H216O entries are 0.0 (non-fractionating)
  static constexpr Scalar aksmc[num_species] = {
    sp(0.0),      // H216O
    sp(0.00528),  // HDO
    sp(0.006),    // H218O
    sp(0.006),    // H217O (assumed same as H218O)
    sp(0.00528)   // HTO (assumed same as HDO)
  };

  // Parameter A for kinetic fractionation
  static constexpr Scalar akrfa[num_species] = {
    sp(0.0),        // H216O
    sp(0.2508e-3),  // HDO
    sp(0.285e-3),   // H218O
    sp(0.285e-3),   // H217O (assumed same as H218O)
    sp(0.2508e-3)   // HTO (assumed same as HDO)
  };

  // Parameter B for kinetic fractionation
  static constexpr Scalar akrfb[num_species] = {
    sp(0.0),        // H216O
    sp(0.7216e-3),  // HDO
    sp(0.82e-3),    // H218O
    sp(0.82e-3),    // H217O (assumed same as H218O)
    sp(0.7216e-3)   // HTO (assumed same as HDO)
  };

  // -----------------------------------------------------------------------
  // Liquid/vapor equilibrium fractionation coefficients
  // -----------------------------------------------------------------------

  // Coefficients for alpha = exp(polynomial in T)
  // Default: Horita & Wesolowski 1994 (Fortran lines 157-163)
  // Alternative: WISO_USE_MAJOUBE_1971 for Majoube 1971a
  //
  // HDO uses: alpha = exp(a*T^3 + b*T^2 + c*T + d + e/T^3)
  // H218O uses: alpha = exp(a/T^3 + b/T^2 + c/T + d)
  // H217O and HTO are derived from these by power laws

#ifdef WISO_USE_MAJOUBE_1971
  // Majoube 1971a (Fortran lines 151-155)
  // Coefficient A
  static constexpr Scalar alpal[num_species] = {
    sp(0.0),       // H216O
    sp(24.844e3),  // HDO
    sp(1.137e3),   // H218O
    sp(0.0),       // H217O (placeholder)
    sp(0.0)        // HTO (placeholder)
  };

  // Coefficient B
  static constexpr Scalar alpbl[num_species] = {
    sp(0.0),      // H216O
    sp(-76.248),  // HDO
    sp(-0.4156),  // H218O
    sp(0.0),      // H217O (placeholder)
    sp(0.0)       // HTO (placeholder)
  };

  // Coefficient C
  static constexpr Scalar alpcl[num_species] = {
    sp(0.0),         // H216O
    sp(52.612e-3),   // HDO
    sp(-2.0667e-3),  // H218O
    sp(0.0),         // H217O (placeholder)
    sp(0.0)          // HTO (placeholder)
  };

  // Coefficient D (not used in Majoube 1971a)
  static constexpr Scalar alpdl[num_species] = {
    sp(0.0), sp(0.0), sp(0.0), sp(0.0), sp(0.0)
  };

  // Coefficient E (not used in Majoube 1971a)
  static constexpr Scalar alpel[num_species] = {
    sp(0.0), sp(0.0), sp(0.0), sp(0.0), sp(0.0)
  };
#else
  // Horita & Wesolowski 1994 (default, Fortran lines 157-163)
  // Coefficient A
  static constexpr Scalar alpal[num_species] = {
    sp(0.0),         // H216O
    sp(1158.8e-12),  // HDO
    sp(0.35041e6),   // H218O
    sp(0.0),         // H217O (derived from H218O)
    sp(0.0)          // HTO (derived from HDO)
  };

  // Coefficient B
  static constexpr Scalar alpbl[num_species] = {
    sp(0.0),        // H216O
    sp(-1620.1e-9), // HDO
    sp(-1.6664e3),  // H218O
    sp(0.0),        // H217O (derived from H218O)
    sp(0.0)         // HTO (derived from HDO)
  };

  // Coefficient C
  static constexpr Scalar alpcl[num_species] = {
    sp(0.0),       // H216O
    sp(794.84e-6), // HDO
    sp(6.7123),    // H218O
    sp(0.0),       // H217O (derived from H218O)
    sp(0.0)        // HTO (derived from HDO)
  };

  // Coefficient D
  static constexpr Scalar alpdl[num_species] = {
    sp(0.0),       // H216O
    sp(-161.04e-3), // HDO
    sp(-7.685e-3), // H218O
    sp(0.0),       // H217O (derived from H218O)
    sp(0.0)        // HTO (derived from HDO)
  };

  // Coefficient E (only used for HDO)
  static constexpr Scalar alpel[num_species] = {
    sp(0.0),       // H216O
    sp(2.9992e6),  // HDO
    sp(0.0),       // H218O
    sp(0.0),       // H217O (derived from H218O)
    sp(0.0)        // HTO (derived from HDO)
  };
#endif

  // -----------------------------------------------------------------------
  // Ice/vapor equilibrium fractionation coefficients
  // -----------------------------------------------------------------------

  // Coefficients for alpha = exp(a/T^2 + b/T + c)
  // Default: Merlivat & Nief 1967 for HDO, Majoube 1971b for H218O
  //          (Fortran lines 172-175)
  // Alternative: WISO_USE_ISOCAM3_ICE for isoCAM3 values
  //              (Fortran lines 165-169)

#ifdef WISO_USE_ISOCAM3_ICE
  // isoCAM3 values (Fortran lines 165-169)
  // Coefficient A
  static constexpr Scalar alpai[num_species] = {
    sp(0.0),     // H216O
    sp(16288.0), // HDO
    sp(0.0),     // H218O
    sp(0.0),     // H217O (derived from H218O)
    sp(0.0)      // HTO (derived from HDO)
  };

  // Coefficient B
  static constexpr Scalar alpbi[num_species] = {
    sp(0.0),    // H216O
    sp(0.0),    // HDO
    sp(11.839), // H218O
    sp(0.0),    // H217O (derived from H218O)
    sp(0.0)     // HTO (derived from HDO)
  };

  // Coefficient C
  static constexpr Scalar alpci[num_species] = {
    sp(0.0),        // H216O
    sp(-9.34e-2),   // HDO
    sp(-28.224e-3), // H218O
    sp(0.0),        // H217O (derived from H218O)
    sp(0.0)         // HTO (derived from HDO)
  };
#else
  // Merlivat & Nief 1967 (HDO), Majoube 1971b (H218O) - default
  // (Fortran lines 172-175)
  // Coefficient A
  static constexpr Scalar alpai[num_species] = {
    sp(0.0),     // H216O
    sp(16289.0), // HDO (Merlivat & Nief 1967)
    sp(0.0),     // H218O
    sp(0.0),     // H217O (derived from H218O)
    sp(0.0)      // HTO (derived from HDO)
  };

  // Coefficient B
  static constexpr Scalar alpbi[num_species] = {
    sp(0.0),    // H216O
    sp(0.0),    // HDO
    sp(11.839), // H218O (Majoube 1971b)
    sp(0.0),    // H217O (derived from H218O)
    sp(0.0)     // HTO (derived from HDO)
  };

  // Coefficient C
  static constexpr Scalar alpci[num_species] = {
    sp(0.0),        // H216O
    sp(-9.45e-2),   // HDO (Merlivat & Nief 1967)
    sp(-28.224e-3), // H218O (Majoube 1971b)
    sp(0.0),        // H217O (derived from H218O)
    sp(0.0)         // HTO (derived from HDO)
  };
#endif

  // -----------------------------------------------------------------------
  // Physical constants for kinetic calculations
  // -----------------------------------------------------------------------

  // Molecular diffusivity of air [m2/s]
  // From water_isotopes.F90 line 215
  static constexpr Scalar difair = sp(2.36e-5);

  // Dynamic viscosity of air [Pa*s or kg/(m*s)]
  // From water_isotopes.F90 line 216
  // Note: about 17°C, 1.73e-5 at STP (Salby)
  static constexpr Scalar muair = sp(1.7e-5);

  // -----------------------------------------------------------------------
  // Safety parameters
  // -----------------------------------------------------------------------

  // Small value for safe division (avoid division by zero)
  // From water_isotopes.F90 line 689
  static constexpr Scalar qtiny = sp(1.e-16);

  // Minimum water threshold for tracer calculations [kg/kg]
  // From eamxx_water_tracers_functions.hpp line 34
  static constexpr Scalar wtrc_qmin = sp(1.e-22);
};

// Convenience alias for Real precision
using WaterIsotopeConstantsReal = WaterIsotopeConstants<Real>;

} // namespace wiso
} // namespace scream

#endif // EAMXX_WATER_ISOTOPES_CONSTANTS_HPP
