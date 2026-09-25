#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_constants.hpp"
#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_fractionation.hpp"

#include "share/core/eamxx_types.hpp"

#include <cmath>

namespace scream {
namespace {

using wiso::WaterIsotopeConstants;
using wiso::WaterIsotopeRuntimeOptions;
using wiso::WaterIsotopeFractionation;

TEST_CASE("runtime_formulation_selection") {
  using Real = scream::Real;

  SECTION("standard_ratio_formulations") {
    // Test 1: Default (Normalized)
    WaterIsotopeRuntimeOptions opts_normalized;
    WaterIsotopeConstants<Real> constants_normalized(opts_normalized);

    // Verify all 1.0
    for (int i = 0; i < 5; ++i) {
      REQUIRE(constants_normalized.ratio_src(static_cast<wiso::WaterIsotopologues>(i)) == Real(1.0));
    }

    // Test 2: Natural abundance
    WaterIsotopeRuntimeOptions opts_natural;
    opts_natural.standard_ratio = wiso::StandardRatioFormulation::NaturalAbundance;
    WaterIsotopeConstants<Real> constants_natural(opts_natural);

    // Verify natural abundance values
    REQUIRE(constants_natural.ratio_src(wiso::WaterIsotopologues::H216O) == Real(1.0));
    REQUIRE(std::abs(constants_natural.ratio_src(wiso::WaterIsotopologues::HDO) - Real(155.76e-6)) < 1e-8);
    REQUIRE(std::abs(constants_natural.ratio_src(wiso::WaterIsotopologues::H218O) - Real(2005.20e-6)) < 1e-6);

    // Verify they differ
    REQUIRE(constants_normalized.ratio_src(wiso::WaterIsotopologues::HDO) != constants_natural.ratio_src(wiso::WaterIsotopologues::HDO));
  }

  SECTION("ocean_enrichment_formulations") {
    // Test 1: Default (None)
    WaterIsotopeRuntimeOptions opts_none;
    WaterIsotopeConstants<Real> constants_none(opts_none);

    // Verify all 1.0
    for (int i = 0; i < 5; ++i) {
      REQUIRE(constants_none.ocean_src(static_cast<wiso::WaterIsotopologues>(i)) == Real(1.0));
    }

    // Test 2: LGM
    WaterIsotopeRuntimeOptions opts_lgm;
    opts_lgm.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    WaterIsotopeConstants<Real> constants_lgm(opts_lgm);

    // Verify LGM values
    REQUIRE(constants_lgm.ocean_src(wiso::WaterIsotopologues::H216O) == Real(1.0));
    REQUIRE(std::abs(constants_lgm.ocean_src(wiso::WaterIsotopologues::HDO) - Real(1.0128)) < 1e-5);
    REQUIRE(std::abs(constants_lgm.ocean_src(wiso::WaterIsotopologues::H218O) - Real(1.0016)) < 1e-5);

    // Verify they differ
    REQUIRE(constants_none.ocean_src(wiso::WaterIsotopologues::HDO) != constants_lgm.ocean_src(wiso::WaterIsotopologues::HDO));
  }

  SECTION("liquid_vapor_fractionation") {
    using wiso::CondensedPhase;
    using wiso::IsoElement;

    // Test 1: Default (Horita & Wesolowski 1994)
    WaterIsotopeRuntimeOptions opts_horita;
    WaterIsotopeConstants<Real> constants_horita(opts_horita);

    // Verify Horita & Wesolowski 1994 hydrogen coefficients. NOTE the per-mil
    // convention: the table holds 10^3 * ln(alpha), so each value is 1000x the
    // corresponding coefficient of ln(alpha).
    const auto& horita_h = constants_horita.alpha_eq_coeffs(CondensedPhase::Liquid,
                                                            IsoElement::Hydrogen);
    REQUIRE(std::abs(horita_h.T3 - Real(1.1588e-6)) < 1e-12);
    REQUIRE(std::abs(horita_h.T_3 - Real(2.9992e9)) < 1e2);

    // Test 2: Majoube 1971
    WaterIsotopeRuntimeOptions opts_majoube;
    opts_majoube.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    WaterIsotopeConstants<Real> constants_majoube(opts_majoube);

    const auto& majoube_h = constants_majoube.alpha_eq_coeffs(CondensedPhase::Liquid,
                                                               IsoElement::Hydrogen);
    REQUIRE(std::abs(majoube_h.T_2 - Real(2.4844e7)) < 1e0);
    // Majoube is a pure function of 1/T: no ascending-power terms.
    REQUIRE(majoube_h.T3 == Real(0.0));
    REQUIRE(majoube_h.T2 == Real(0.0));
    REQUIRE(majoube_h.T1 == Real(0.0));

    // Verify they differ significantly
    REQUIRE(std::abs(horita_h.T_3 - majoube_h.T_3) > 1e6);

    // Fitted ranges differ: Horita extends to 364 C, Majoube stops at 100 C.
    REQUIRE(constants_horita.tbounds(CondensedPhase::Liquid, IsoElement::Hydrogen).Tmax >
            constants_majoube.tbounds(CondensedPhase::Liquid, IsoElement::Hydrogen).Tmax);
  }

  SECTION("ice_vapor_fractionation") {
    using wiso::CondensedPhase;
    using wiso::IsoElement;

    // Test 1: Default (Merlivat & Nief 1967 for HDO, Majoube 1971 for H218O)
    WaterIsotopeRuntimeOptions opts_merlivat;
    WaterIsotopeConstants<Real> constants_merlivat(opts_merlivat);

    const auto& merlivat_h = constants_merlivat.alpha_eq_coeffs(CondensedPhase::Ice,
                                                                 IsoElement::Hydrogen);
    REQUIRE(std::abs(merlivat_h.T_2 - Real(1.6289e7)) < 1e0);
    REQUIRE(std::abs(merlivat_h.T0 - Real(-9.45e1)) < 1e-1);

    // The oxygen row comes from a different study, with a different fitted
    // range: 1/T only, and starting at 239.75 K rather than 233.15 K.
    const auto& majoube_o = constants_merlivat.alpha_eq_coeffs(CondensedPhase::Ice,
                                                               IsoElement::Oxygen);
    REQUIRE(std::abs(majoube_o.T_1 - Real(1.1839e4)) < 1e0);
    REQUIRE(majoube_o.T_2 == Real(0.0));
    REQUIRE(std::abs(constants_merlivat.tbounds(CondensedPhase::Ice, IsoElement::Oxygen).Tmin
                     - Real(239.75)) < 1e-2);

    // Test 2: isoCAM3
    WaterIsotopeRuntimeOptions opts_isocam3;
    opts_isocam3.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;
    WaterIsotopeConstants<Real> constants_isocam3(opts_isocam3);

    const auto& isocam3_h = constants_isocam3.alpha_eq_coeffs(CondensedPhase::Ice,
                                                               IsoElement::Hydrogen);
    REQUIRE(std::abs(isocam3_h.T_2 - Real(1.6288e7)) < 1e0);
    REQUIRE(std::abs(isocam3_h.T0 - Real(-9.34e1)) < 1e-1);

    // Verify they differ (slightly, but measurably)
    REQUIRE(merlivat_h.T_2 != isocam3_h.T_2);
    REQUIRE(merlivat_h.T0 != isocam3_h.T0);

    // isoCAM3 exists to extrapolate colder than the original fits.
    REQUIRE(constants_isocam3.tbounds(CondensedPhase::Ice, IsoElement::Hydrogen).Tmin <
            constants_merlivat.tbounds(CondensedPhase::Ice, IsoElement::Hydrogen).Tmin);
  }

  SECTION("combined_formulations") {
    // Test using all alternative formulations together
    WaterIsotopeRuntimeOptions opts_alt;
    opts_alt.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    opts_alt.standard_ratio = wiso::StandardRatioFormulation::NaturalAbundance;
    opts_alt.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    opts_alt.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;

    WaterIsotopeConstants<Real> constants_alt(opts_alt);
    WaterIsotopeConstants<Real> constants_default;  // All defaults

    // Verify each setting was applied correctly
    REQUIRE(constants_alt.ratio_src(wiso::WaterIsotopologues::HDO) != constants_default.ratio_src(wiso::WaterIsotopologues::HDO));
    REQUIRE(constants_alt.ocean_src(wiso::WaterIsotopologues::HDO) != constants_default.ocean_src(wiso::WaterIsotopologues::HDO));
    REQUIRE(constants_alt.alpha_eq_coeffs(wiso::CondensedPhase::Liquid, wiso::IsoElement::Hydrogen).T_3 !=
            constants_default.alpha_eq_coeffs(wiso::CondensedPhase::Liquid, wiso::IsoElement::Hydrogen).T_3);
    REQUIRE(constants_alt.alpha_eq_coeffs(wiso::CondensedPhase::Ice, wiso::IsoElement::Hydrogen).T_2 !=
            constants_default.alpha_eq_coeffs(wiso::CondensedPhase::Ice, wiso::IsoElement::Hydrogen).T_2);
  }

  SECTION("fractionation_with_runtime_constants") {
    // Test that fractionation functions produce different results with different formulations
    const Real temp = Real(273.15);  // 0°C

    // Horita & Wesolowski 1994
    WaterIsotopeRuntimeOptions opts_horita;
    WaterIsotopeConstants<Real> constants_horita(opts_horita);
    Real alpha_lv_horita = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, constants_horita);

    // Majoube 1971
    WaterIsotopeRuntimeOptions opts_majoube;
    opts_majoube.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    WaterIsotopeConstants<Real> constants_majoube(opts_majoube);
    Real alpha_lv_majoube = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, constants_majoube);

    // Verify they produce different fractionation factors
    REQUIRE(alpha_lv_horita != alpha_lv_majoube);
    REQUIRE(alpha_lv_horita > Real(1.0));  // Should enrich heavy isotope
    REQUIRE(alpha_lv_majoube > Real(1.0));

    // Relative difference should be measurable
    Real rel_diff = std::abs(alpha_lv_horita - alpha_lv_majoube) / alpha_lv_horita;
    REQUIRE(rel_diff > 1e-6);  // At least 0.0001% difference
  }

  SECTION("phase_wrappers_match_generic_form") {
    // The phase-specific spellings must be exactly the generic evaluator with
    // the phase pinned -- they add no logic, so they must not add any drift.
    const Real temp = Real(263.15);
    WaterIsotopeConstants<Real> constants;

    REQUIRE(WaterIsotopeFractionation::alpha_liquid_vapor(
              temp, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, constants) ==
            WaterIsotopeFractionation::alpha_equilibrium(
              temp, wiso::WaterIsotopologues::HDO, wiso::CondensedPhase::Liquid,
              wiso::CondensedOverVapor, constants));

    REQUIRE(WaterIsotopeFractionation::alpha_ice_vapor(
              temp, wiso::WaterIsotopologues::H218O, wiso::CondensedOverVapor, constants) ==
            WaterIsotopeFractionation::alpha_equilibrium(
              temp, wiso::WaterIsotopologues::H218O, wiso::CondensedPhase::Ice,
              wiso::CondensedOverVapor, constants));
  }
}

} // anonymous namespace
} // namespace scream
