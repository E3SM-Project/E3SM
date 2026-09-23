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

  SECTION("diffusivity_formulations") {
    // Test 1: Default (Merlivat 1978)
    WaterIsotopeRuntimeOptions opts_default;
    WaterIsotopeConstants<Real> constants_default(opts_default);

    // Verify Merlivat 1978 values
    REQUIRE(constants_default.diff_src(wiso::H216O) == Real(1.0));
    REQUIRE(std::abs(constants_default.diff_src(wiso::HDO) - Real(0.9757)) < 1e-5);
    REQUIRE(std::abs(constants_default.diff_src(wiso::H218O) - Real(0.9727)) < 1e-5);

    // Test 2: Cappa 2003
    WaterIsotopeRuntimeOptions opts_cappa;
    opts_cappa.diffusivity = wiso::DiffusivityFormulation::Cappa2003;
    WaterIsotopeConstants<Real> constants_cappa(opts_cappa);

    // Verify Cappa 2003 values
    REQUIRE(constants_cappa.diff_src(wiso::H216O) == Real(1.0));
    REQUIRE(std::abs(constants_cappa.diff_src(wiso::HDO) - Real(0.9839)) < 1e-5);
    REQUIRE(std::abs(constants_cappa.diff_src(wiso::H218O) - Real(0.9691)) < 1e-5);

    // Verify they differ
    REQUIRE(constants_default.diff_src(wiso::HDO) != constants_cappa.diff_src(wiso::HDO));
    REQUIRE(constants_default.diff_src(wiso::H218O) != constants_cappa.diff_src(wiso::H218O));
  }

  SECTION("standard_ratio_formulations") {
    // Test 1: Default (Normalized)
    WaterIsotopeRuntimeOptions opts_normalized;
    WaterIsotopeConstants<Real> constants_normalized(opts_normalized);

    // Verify all 1.0
    for (int i = 0; i < 5; ++i) {
      REQUIRE(constants_normalized.ratio_src(i) == Real(1.0));
    }

    // Test 2: Natural abundance
    WaterIsotopeRuntimeOptions opts_natural;
    opts_natural.standard_ratio = wiso::StandardRatioFormulation::NaturalAbundance;
    WaterIsotopeConstants<Real> constants_natural(opts_natural);

    // Verify natural abundance values
    REQUIRE(constants_natural.ratio_src(wiso::H216O) == Real(1.0));
    REQUIRE(std::abs(constants_natural.ratio_src(wiso::HDO) - Real(155.76e-6)) < 1e-8);
    REQUIRE(std::abs(constants_natural.ratio_src(wiso::H218O) - Real(2005.20e-6)) < 1e-6);

    // Verify they differ
    REQUIRE(constants_normalized.ratio_src(wiso::HDO) != constants_natural.ratio_src(wiso::HDO));
  }

  SECTION("ocean_enrichment_formulations") {
    // Test 1: Default (None)
    WaterIsotopeRuntimeOptions opts_none;
    WaterIsotopeConstants<Real> constants_none(opts_none);

    // Verify all 1.0
    for (int i = 0; i < 5; ++i) {
      REQUIRE(constants_none.ocean_src(i) == Real(1.0));
    }

    // Test 2: LGM
    WaterIsotopeRuntimeOptions opts_lgm;
    opts_lgm.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    WaterIsotopeConstants<Real> constants_lgm(opts_lgm);

    // Verify LGM values
    REQUIRE(constants_lgm.ocean_src(wiso::H216O) == Real(1.0));
    REQUIRE(std::abs(constants_lgm.ocean_src(wiso::HDO) - Real(1.0128)) < 1e-5);
    REQUIRE(std::abs(constants_lgm.ocean_src(wiso::H218O) - Real(1.0016)) < 1e-5);

    // Verify they differ
    REQUIRE(constants_none.ocean_src(wiso::HDO) != constants_lgm.ocean_src(wiso::HDO));
  }

  SECTION("liquid_vapor_fractionation") {
    // Test 1: Default (Horita & Wesolowski 1994)
    WaterIsotopeRuntimeOptions opts_horita;
    WaterIsotopeConstants<Real> constants_horita(opts_horita);

    // Verify Horita & Wesolowski 1994 values for HDO coefficient A
    REQUIRE(std::abs(constants_horita.AlphaLiqVap_CoefA(wiso::HDO) - Real(1158.8e-12)) < 1e-16);
    REQUIRE(std::abs(constants_horita.AlphaLiqVap_CoefE(wiso::HDO) - Real(2.9992e6)) < 1e-1);

    // Test 2: Majoube 1971
    WaterIsotopeRuntimeOptions opts_majoube;
    opts_majoube.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    WaterIsotopeConstants<Real> constants_majoube(opts_majoube);

    // Verify Majoube 1971a values for HDO coefficient A
    REQUIRE(std::abs(constants_majoube.AlphaLiqVap_CoefA(wiso::HDO) - Real(24.844e3)) < 1e-1);
    // Majoube doesn't use E coefficient
    REQUIRE(constants_majoube.AlphaLiqVap_CoefE(wiso::HDO) == Real(0.0));

    // Verify they differ significantly
    REQUIRE(std::abs(constants_horita.AlphaLiqVap_CoefA(wiso::HDO) - constants_majoube.AlphaLiqVap_CoefA(wiso::HDO)) > 1e3);
  }

  SECTION("ice_vapor_fractionation") {
    // Test 1: Default (Merlivat & Nief 1967)
    WaterIsotopeRuntimeOptions opts_merlivat;
    WaterIsotopeConstants<Real> constants_merlivat(opts_merlivat);

    // Verify Merlivat & Nief 1967 values for HDO coefficient A
    REQUIRE(std::abs(constants_merlivat.AlphaIceVap_CoefA(wiso::HDO) - Real(16289.0)) < 1e-1);
    REQUIRE(std::abs(constants_merlivat.AlphaIceVap_CoefC(wiso::HDO) - Real(-9.45e-2)) < 1e-4);

    // Test 2: isoCAM3
    WaterIsotopeRuntimeOptions opts_isocam3;
    opts_isocam3.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;
    WaterIsotopeConstants<Real> constants_isocam3(opts_isocam3);

    // Verify isoCAM3 values for HDO coefficient A
    REQUIRE(std::abs(constants_isocam3.AlphaIceVap_CoefA(wiso::HDO) - Real(16288.0)) < 1e-1);
    REQUIRE(std::abs(constants_isocam3.AlphaIceVap_CoefC(wiso::HDO) - Real(-9.34e-2)) < 1e-4);

    // Verify they differ (slightly, but measurably)
    REQUIRE(constants_merlivat.AlphaIceVap_CoefA(wiso::HDO) != constants_isocam3.AlphaIceVap_CoefA(wiso::HDO));
    REQUIRE(constants_merlivat.AlphaIceVap_CoefC(wiso::HDO) != constants_isocam3.AlphaIceVap_CoefC(wiso::HDO));
  }

  SECTION("combined_formulations") {
    // Test using all alternative formulations together
    WaterIsotopeRuntimeOptions opts_alt;
    opts_alt.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    opts_alt.diffusivity = wiso::DiffusivityFormulation::Cappa2003;
    opts_alt.standard_ratio = wiso::StandardRatioFormulation::NaturalAbundance;
    opts_alt.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    opts_alt.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;

    WaterIsotopeConstants<Real> constants_alt(opts_alt);
    WaterIsotopeConstants<Real> constants_default;  // All defaults

    // Verify each setting was applied correctly
    REQUIRE(constants_alt.diff_src(wiso::HDO) != constants_default.diff_src(wiso::HDO));
    REQUIRE(constants_alt.ratio_src(wiso::HDO) != constants_default.ratio_src(wiso::HDO));
    REQUIRE(constants_alt.ocean_src(wiso::HDO) != constants_default.ocean_src(wiso::HDO));
    REQUIRE(constants_alt.AlphaLiqVap_CoefA(wiso::HDO) != constants_default.AlphaLiqVap_CoefA(wiso::HDO));
    REQUIRE(constants_alt.AlphaIceVap_CoefA(wiso::HDO) != constants_default.AlphaIceVap_CoefA(wiso::HDO));
  }

  SECTION("fractionation_with_runtime_constants") {
    // Test that fractionation functions produce different results with different formulations
    const Real temp = Real(273.15);  // 0°C

    // Horita & Wesolowski 1994
    WaterIsotopeRuntimeOptions opts_horita;
    WaterIsotopeConstants<Real> constants_horita(opts_horita);
    Real alpha_lv_horita = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::HDO, wiso::CondensedOverVapor, constants_horita);

    // Majoube 1971
    WaterIsotopeRuntimeOptions opts_majoube;
    opts_majoube.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    WaterIsotopeConstants<Real> constants_majoube(opts_majoube);
    Real alpha_lv_majoube = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::HDO, wiso::CondensedOverVapor, constants_majoube);

    // Verify they produce different fractionation factors
    REQUIRE(alpha_lv_horita != alpha_lv_majoube);
    REQUIRE(alpha_lv_horita > Real(1.0));  // Should enrich heavy isotope
    REQUIRE(alpha_lv_majoube > Real(1.0));

    // Relative difference should be measurable
    Real rel_diff = std::abs(alpha_lv_horita - alpha_lv_majoube) / alpha_lv_horita;
    REQUIRE(rel_diff > 1e-6);  // At least 0.0001% difference
  }

  SECTION("backward_compatibility") {
    // Test that 3-argument form still works and matches default formulation
    const Real temp = Real(273.15);

    // Call 3-argument form (backward compatible)
    Real alpha_3arg = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::HDO, wiso::CondensedOverVapor);

    // Call 4-argument form with default constants
    WaterIsotopeConstants<Real> constants_default;
    Real alpha_4arg = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::HDO, wiso::CondensedOverVapor, constants_default);

    // They should be identical
    REQUIRE(alpha_3arg == alpha_4arg);
  }
}

} // anonymous namespace
} // namespace scream
