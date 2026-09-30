#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_parameters.hpp"
#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_fractionation.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_pack.hpp>

#include <cmath>

namespace scream {
namespace {

using wiso::WaterIsotopeParameters;
using wiso::WaterIsotopeRuntimeOptions;
using wiso::WaterIsotopeFractionation;

TEST_CASE("runtime_formulation_selection") {

  SECTION("standard_ratio_formulations") {
    // Test 1: Default (Normalized)
    WaterIsotopeRuntimeOptions opts_normalized;
    WaterIsotopeParameters<Real> iso_params_normalized(opts_normalized);

    // Verify all 1.0
    for (int i = 0; i < 5; ++i) {
      REQUIRE(iso_params_normalized.std_iso_ratio(static_cast<wiso::WaterIsotopologues>(i)) == Real(1.0));
    }

    // Test 2: Natural abundance
    WaterIsotopeRuntimeOptions opts_natural;
    opts_natural.standard_ratio = wiso::StandardRatioFormulation::NaturalAbundance;
    WaterIsotopeParameters<Real> iso_params_natural(opts_natural);

    // Verify natural abundance values
    REQUIRE(iso_params_natural.std_iso_ratio(wiso::WaterIsotopologues::H216O) == Real(0.9976));
    REQUIRE(std::abs(iso_params_natural.std_iso_ratio(wiso::WaterIsotopologues::HDO) - Real(155.76e-6)) < 1e-8);
    REQUIRE(std::abs(iso_params_natural.std_iso_ratio(wiso::WaterIsotopologues::H218O) - Real(2005.20e-6)) < 1e-6);

    // Verify they differ
    REQUIRE(iso_params_normalized.std_iso_ratio(wiso::WaterIsotopologues::HDO) != iso_params_natural.std_iso_ratio(wiso::WaterIsotopologues::HDO));
  }

  SECTION("ocean_enrichment_formulations") {
    // Test 1: Default (None)
    WaterIsotopeRuntimeOptions opts_none;
    WaterIsotopeParameters<Real> iso_params_none(opts_none);

    // Verify all 1.0
    for (int i = 0; i < 5; ++i) {
      REQUIRE(iso_params_none.mean_ocean_enrichment(static_cast<wiso::WaterIsotopologues>(i)) == Real(1.0));
    }

    // Test 2: LGM
    WaterIsotopeRuntimeOptions opts_lgm;
    opts_lgm.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    WaterIsotopeParameters<Real> iso_params_lgm(opts_lgm);

    // Verify LGM values
    REQUIRE(iso_params_lgm.mean_ocean_enrichment(wiso::WaterIsotopologues::H216O) == Real(1.0));
    REQUIRE(std::abs(iso_params_lgm.mean_ocean_enrichment(wiso::WaterIsotopologues::HDO) - Real(1.0128)) < 1e-5);
    REQUIRE(std::abs(iso_params_lgm.mean_ocean_enrichment(wiso::WaterIsotopologues::H218O) - Real(1.0016)) < 1e-5);

    // Verify they differ
    REQUIRE(iso_params_none.mean_ocean_enrichment(wiso::WaterIsotopologues::HDO) != iso_params_lgm.mean_ocean_enrichment(wiso::WaterIsotopologues::HDO));
  }

  SECTION("liquid_vapor_fractionation") {
    using wiso::CondensedPhase;
    using wiso::IsoElement;

    // Test 1: Default (Horita & Wesolowski 1994)
    WaterIsotopeRuntimeOptions opts_horita;
    WaterIsotopeParameters<Real> iso_params_horita(opts_horita);

    // Verify Horita & Wesolowski 1994 hydrogen coefficients. NOTE the per-mil
    // convention: the table holds 10^3 * ln(alpha), so each value is 1000x the
    // corresponding coefficient of ln(alpha).
    // RPF note - do these tests add much? Consider rethinking them.
    const auto& horita_h = iso_params_horita.alpha_eq_coeffs(CondensedPhase::Liquid,
                                                            IsoElement::Hydrogen);
    REQUIRE(std::abs(horita_h.T3 - Real(1.1588e-6)) < 1e-12);
    REQUIRE(std::abs(horita_h.T_3 - Real(2.9992e9)) < 1e2);

    // Test 2: Majoube 1971
    WaterIsotopeRuntimeOptions opts_majoube;
    opts_majoube.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    WaterIsotopeParameters<Real> iso_params_majoube(opts_majoube);

    const auto& majoube_h = iso_params_majoube.alpha_eq_coeffs(CondensedPhase::Liquid,
                                                               IsoElement::Hydrogen);
    REQUIRE(std::abs(majoube_h.T_2 - Real(2.4844e7)) < 1e0);
    // Majoube is a pure function of 1/T: no ascending-power terms.
    REQUIRE(majoube_h.T3 == Real(0.0));
    REQUIRE(majoube_h.T2 == Real(0.0));
    REQUIRE(majoube_h.T1 == Real(0.0));

    // Verify they differ significantly
    REQUIRE(std::abs(horita_h.T_3 - majoube_h.T_3) > 1e6);

    // Fitted ranges differ: Horita extends to 364 C, Majoube stops at 100 C.
    REQUIRE(iso_params_horita.tbounds(CondensedPhase::Liquid, IsoElement::Hydrogen).Tmax >
            iso_params_majoube.tbounds(CondensedPhase::Liquid, IsoElement::Hydrogen).Tmax);
  }

  SECTION("ice_vapor_fractionation") {
    using wiso::CondensedPhase;
    using wiso::IsoElement;

    // Test 1: Default (Merlivat & Nief 1967 for HDO, Majoube 1971 for H218O)
    WaterIsotopeRuntimeOptions opts_merlivat;
    WaterIsotopeParameters<Real> iso_params_merlivat(opts_merlivat);

    const auto& merlivat_h = iso_params_merlivat.alpha_eq_coeffs(CondensedPhase::Ice,
                                                                 IsoElement::Hydrogen);
    REQUIRE(std::abs(merlivat_h.T_2 - Real(1.6289e7)) < 1e0);
    REQUIRE(std::abs(merlivat_h.T0 - Real(-9.45e1)) < 1e-1);

    // The oxygen row comes from a different study, with a different fitted
    // range: 1/T only, and starting at 239.75 K rather than 233.15 K.
    const auto& majoube_o = iso_params_merlivat.alpha_eq_coeffs(CondensedPhase::Ice,
                                                               IsoElement::Oxygen);
    REQUIRE(std::abs(majoube_o.T_1 - Real(1.1839e4)) < 1e0);
    REQUIRE(majoube_o.T_2 == Real(0.0));
    REQUIRE(std::abs(iso_params_merlivat.tbounds(CondensedPhase::Ice, IsoElement::Oxygen).Tmin
                     - Real(239.75)) < 1e-2);

    // Test 2: isoCAM3
    WaterIsotopeRuntimeOptions opts_isocam3;
    opts_isocam3.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;
    WaterIsotopeParameters<Real> iso_params_isocam3(opts_isocam3);

    const auto& isocam3_h = iso_params_isocam3.alpha_eq_coeffs(CondensedPhase::Ice,
                                                               IsoElement::Hydrogen);
    REQUIRE(std::abs(isocam3_h.T_2 - Real(1.6288e7)) < 1e0);
    REQUIRE(std::abs(isocam3_h.T0 - Real(-9.34e1)) < 1e-1);

    // Verify they differ (slightly, but measurably)
    REQUIRE(merlivat_h.T_2 != isocam3_h.T_2);
    REQUIRE(merlivat_h.T0 != isocam3_h.T0);

    // isoCAM3 exists to extrapolate colder than the original fits.
    REQUIRE(iso_params_isocam3.tbounds(CondensedPhase::Ice, IsoElement::Hydrogen).Tmin <
            iso_params_merlivat.tbounds(CondensedPhase::Ice, IsoElement::Hydrogen).Tmin);
  }

  SECTION("combined_formulations") {
    // Test using all alternative formulations together
    WaterIsotopeRuntimeOptions opts_alt;
    opts_alt.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    opts_alt.standard_ratio = wiso::StandardRatioFormulation::NaturalAbundance;
    opts_alt.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    opts_alt.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;

    WaterIsotopeParameters<Real> iso_params_alt(opts_alt);
    WaterIsotopeParameters<Real> iso_params_default;  // All defaults

    // Verify each setting was applied correctly
    REQUIRE(iso_params_alt.std_iso_ratio(wiso::WaterIsotopologues::HDO) != iso_params_default.std_iso_ratio(wiso::WaterIsotopologues::HDO));
    REQUIRE(iso_params_alt.mean_ocean_enrichment(wiso::WaterIsotopologues::HDO) != iso_params_default.mean_ocean_enrichment(wiso::WaterIsotopologues::HDO));
    REQUIRE(iso_params_alt.alpha_eq_coeffs(wiso::CondensedPhase::Liquid, wiso::IsoElement::Hydrogen).T_3 !=
            iso_params_default.alpha_eq_coeffs(wiso::CondensedPhase::Liquid, wiso::IsoElement::Hydrogen).T_3);
    REQUIRE(iso_params_alt.alpha_eq_coeffs(wiso::CondensedPhase::Ice, wiso::IsoElement::Hydrogen).T_2 !=
            iso_params_default.alpha_eq_coeffs(wiso::CondensedPhase::Ice, wiso::IsoElement::Hydrogen).T_2);
  }

  SECTION("fractionation_with_runtime_constants") {
    // Test that fractionation functions produce different results with different formulations
    using Pack1 = ekat::Pack<Real,1>;
    const Pack1 temp = Pack1(273.15);  // 0°C

    // Horita & Wesolowski 1994
    WaterIsotopeRuntimeOptions opts_horita;
    WaterIsotopeParameters<Real> iso_params_horita(opts_horita);
    Real alpha_lv_horita = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, iso_params_horita)[0];

    // Majoube 1971
    WaterIsotopeRuntimeOptions opts_majoube;
    opts_majoube.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    WaterIsotopeParameters<Real> iso_params_majoube(opts_majoube);
    Real alpha_lv_majoube = WaterIsotopeFractionation::alpha_liquid_vapor(
      temp, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, iso_params_majoube)[0];

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
    using Pack1 = ekat::Pack<Real,1>;
    const Pack1 temp = Pack1(263.15);
    WaterIsotopeParameters<Real> iso_params;

    REQUIRE((WaterIsotopeFractionation::alpha_liquid_vapor(
              temp, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, iso_params) ==
            WaterIsotopeFractionation::alpha_equilibrium(
              temp, wiso::WaterIsotopologues::HDO, wiso::CondensedPhase::Liquid,
              wiso::CondensedOverVapor, iso_params)).all());

    REQUIRE((WaterIsotopeFractionation::alpha_ice_vapor(
              temp, wiso::WaterIsotopologues::H218O, wiso::CondensedOverVapor, iso_params) ==
            WaterIsotopeFractionation::alpha_equilibrium(
              temp, wiso::WaterIsotopologues::H218O, wiso::CondensedPhase::Ice,
              wiso::CondensedOverVapor, iso_params)).all());
  }
}

} // anonymous namespace
} // namespace scream
