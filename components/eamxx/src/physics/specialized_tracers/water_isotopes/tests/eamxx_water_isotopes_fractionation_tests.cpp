#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_fractionation.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_pack.hpp>
#include <ekat_view_utils.hpp>
#include <ekat_fpe.hpp>

#include <cfenv>
#include <cmath>

namespace scream {
namespace {

template <typename ScalarT>
bool relative_approx(const ScalarT& computed, 
                     double expected, 
                     typename ekat::ScalarTraits<ScalarT>::scalar_type tol)
{                    
  const auto rel_err = ekat::abs((computed - expected) / expected);
  return (rel_err < tol).all();  // Check ALL lanes, not just [0]
} 

/* Expected alpha values, transcribed from fractionation_factors.xlsx (columns
   Q-X, the "Temperature checks (Kelvin)" block). Each entry below is one
   coefficient set as the spreadsheet tabulates it: the reference row it came
   from, the temperatures at which that row is inside its own fitted range, and
   the resulting R_condensed/R_vapor (>= 1).

   Holding the spreadsheet's own numbers -- rather than re-transcribing the
   polynomial into C++ a second time -- keeps the check external to the code
   under test. The spreadsheet evaluates the coefficients in expanded form
   (sum of A*T^k) while the implementation uses Horner, so the tolerance need
   only absorb two evaluation orders of the same polynomial, not any physics
   uncertainty.

   The tables are indexed by the formulation enum and sized by its
   FormulationCount, so adding a formulation extends the sweep automatically:
   the test iterates every enumerator rather than naming them. A new enumerator
   with no row here fails `expected_is_populated` with a pointer to this
   comment, rather than silently going unexercised.

   Temperatures are whole Kelvin because the spreadsheet's range test is
   C_column + 273, not + 273.15. The T = 273 columns therefore land 0.15 K
   below the implementation's liquid-phase Tmin of 273.15 K and will trip the
   debug-only "extrapolating" warning. That is a 0.15 K disagreement about
   where 0 C is, not a failing check.

   To add or regenerate a row: read columns Q-X of the spreadsheet row named in
   `reference`, skipping any cell that reads "WARN" (out of range). The Hydrogen
   and Oxygen rows of a formulation must contribute the same temperature list,
   so keep only the columns in range for both. */
constexpr int MAX_TEMPS = 8;

struct AlphaTable {
  const char* reference = nullptr;  // spreadsheet row this came from
  int n = 0;                        // number of in-range temperatures
  double T[MAX_TEMPS] = {};         // [K]
  double hdo[MAX_TEMPS] = {};       // alpha for HDO      (Hydrogen coeff row)
  double o18[MAX_TEMPS] = {};       // alpha for H2(18)O  (Oxygen coeff row)
};

// Indexed by LiquidVaporFractionation.
const AlphaTable
liq_expected[etoi(wiso::LiquidVaporFractionation::FormulationCount)] = {
  // Horita & Wesolowski (1994). Spreadsheet rows 12 (D) and 7 (18O).
  { "Horita & Wesolowski 1994, L-V", 5,
    {273.0, 283.0, 293.0, 303.0, 400.0},
    {1.1120343706807947, 1.0971744958332408, 1.0845303130660782,
     1.073696101136939,  1.0188817475908205},
    {1.0118347780075725, 1.0107441038943008, 1.0097913949759019,
     1.0089533859591724, 1.004164554003981} },
  // Majoube (1971). Spreadsheet rows 14 (D) and 8 (18O).
  { "Majoube 1971, L-V", 4,
    {273.0, 283.0, 293.0, 303.0},
    {1.1125581992610687, 1.0978890989567374, 1.0852080906796728, 1.0741973311576687},
    {1.0117350842524433, 1.0107184907529774, 1.0098068293810731, 1.0089862254523621} }
};

// Indexed by IceVaporFractionation.
const AlphaTable
ice_expected[etoi(wiso::IceVaporFractionation::FormulationCount)] = {
  /* Default: the two elements come from two different studies with different
     fitted ranges (Merlivat & Nief 1971 for D, Majoube 1971 for 18O), so this
     row keeps only the temperatures in range for both.
     Spreadsheet rows 17 (D) and 9 (18O). */
  { "Merlivat & Nief 1971 [D] + Majoube 1971 [18O], I-V", 3,
    {248.0, 263.0, 273.0},
    {1.1857133338327439, 1.1514196623562765, 1.1320829093330054},
    {1.0197055439585383, 1.016932973832859,  1.01525752585485} },
  // isoCAM3. Spreadsheet rows 18 (D) and 19 (18O).
  { "isoCAM3, I-V", 4,
    {223.0, 248.0, 263.0, 273.0},
    {1.2638154003170829, 1.1869990364206471, 1.1526702561788347, 1.1333136792477154},
    {1.0251774156985261, 1.0197055439585383, 1.016932973832859,  1.01525752585485} }
};

/* An enumerator added to LiquidVaporFractionation/IceVaporFractionation without
   a matching row above leaves a default-constructed entry. Catch that here so
   the gap reads as "add the expected values" rather than as a vacuous pass. */
bool expected_is_populated(const AlphaTable& t)
{
  return t.reference != nullptr && t.n > 0 && t.n <= MAX_TEMPS;
}

template <typename ScalarT>
void run_sweep(
  const char* phase_name,  // "liquid-vapor" or "ice-vapor"
  const AlphaTable& expected,  // Tabulated reference values
  typename ekat::ScalarTraits<ScalarT>::scalar_type tol,  // Tolerance
  const wiso::WaterIsotopeConstants<typename ekat::ScalarTraits<ScalarT>::scalar_type>& constants)
{
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
  using WIF = wiso::WaterIsotopeFractionation;

  // Choose the alpha function based on phase
  const bool use_liquid_vapor = (std::string(phase_name) == "liquid-vapor");
  auto alpha_fn = [use_liquid_vapor](const ScalarT& t, wiso::WaterIsotopologues species,
                                      wiso::WisoAlphaDir dir,
                                      const wiso::WaterIsotopeConstants<RealT>& constants) -> ScalarT {
    return use_liquid_vapor ? WIF::alpha_liquid_vapor<ScalarT>(t, species, dir, constants)
                            : WIF::alpha_ice_vapor<ScalarT>(t, species, dir, constants);
  };

  INFO("reference: " << expected.reference);

  double prev_hdo = 1e30, prev_o18 = 1e30;
  for (int i = 0; i < expected.n; ++i) {
    const double T_i    = expected.T[i];
    const double ref_hdo = expected.hdo[i];
    const double ref_o18 = expected.o18[i];
    INFO("T = " << T_i << " K");

    const ScalarT t(T_i);

    const ScalarT a_hdo = alpha_fn(t, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor,
constants);
    const ScalarT a_o18 = alpha_fn(t, wiso::WaterIsotopologues::H218O, wiso::CondensedOverVapor,
constants);

    REQUIRE( relative_approx(a_hdo, ref_hdo, tol) );
    REQUIRE( relative_approx(a_o18, ref_o18, tol) );

    // H216O check
    const ScalarT a_16 = alpha_fn(t, wiso::WaterIsotopologues::H216O, wiso::CondensedOverVapor,
constants);
    REQUIRE( (a_16 == ScalarT(1)).all() );
    
    // Direction checks
    const ScalarT a_hdo_inv = alpha_fn(t, wiso::WaterIsotopologues::HDO, wiso::VaporOverCondensed,
constants);
    const ScalarT a_o18_inv = alpha_fn(t, wiso::WaterIsotopologues::H218O, wiso::VaporOverCondensed,
constants);
    REQUIRE( relative_approx(a_hdo_inv, 1.0/ref_hdo, tol) );
    REQUIRE( relative_approx(a_o18_inv, 1.0/ref_o18, tol) );
    
    // Power law checks for H217O and HTO
    const ScalarT a_17 = alpha_fn(t, wiso::WaterIsotopologues::H217O, wiso::CondensedOverVapor,
constants);
    const ScalarT a_ht = alpha_fn(t, wiso::WaterIsotopologues::HTO, wiso::CondensedOverVapor,
constants);
    REQUIRE( relative_approx(a_17, std::pow(ref_o18, 0.529), tol) );
    REQUIRE( relative_approx(a_ht, std::pow(ref_hdo, 2.0), tol) );

    // Monotonicity and >= 1 checks
    REQUIRE( a_hdo[0] >= RealT(1) );
    REQUIRE( a_o18[0] >= RealT(1) );
    REQUIRE( a_hdo[0] < prev_hdo );
    REQUIRE( a_o18[0] < prev_o18 );
    prev_hdo = a_hdo[0];
    prev_o18 = a_o18[0];
    }
}

// Index of a given temperature within a table, so the device check below can
// name the temperature it wants instead of carrying a magic index that would
// silently shift if a spreadsheet column moved in or out of range. Returns -1
// if absent; the caller asserts.
int index_of_T(const AlphaTable& tbl, double T)
{
  for (int i = 0; i < tbl.n; ++i) {
    if (tbl.T[i] == T) return i;
  }
  return -1;
}

// Exercise the KOKKOS_INLINE_FUNCTION on the device to confirm it is
// device-callable and gives the same result as the host path.
void run_on_device()
{
  using WIF = wiso::WaterIsotopeFractionation;
  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using view_1d = typename KT::template view_1d<Real>;

  // The kernel below default-constructs its constants, so check against the
  // default formulations' rows. 273 K is tabulated for both, so one temperature
  // covers both device calls.
  const wiso::WaterIsotopeRuntimeOptions defaults;
  const AlphaTable& liq = liq_expected[etoi(defaults.liquid_vapor)];
  const AlphaTable& ice = ice_expected[etoi(defaults.ice_vapor)];

  const double T_chk = 273.0;
  const int i_liq = index_of_T(liq, T_chk);
  const int i_ice = index_of_T(ice, T_chk);
  REQUIRE( i_liq >= 0 );
  REQUIRE( i_ice >= 0 );

  view_1d out("wiso_alpha_device", 2);
  // Constructed inside the kernel to confirm the resolving constructor itself is
  // device-callable, not just the evaluator.
  Kokkos::parallel_for("wiso_frac_device", 1, KOKKOS_LAMBDA(const int /*i*/) {
    wiso::WaterIsotopeConstants<Real> constants;
    out(0) = WIF::alpha_liquid_vapor(Real(T_chk), wiso::WaterIsotopologues::HDO,   wiso::CondensedOverVapor, constants);
    out(1) = WIF::alpha_ice_vapor   (Real(T_chk), wiso::WaterIsotopologues::H218O, wiso::CondensedOverVapor, constants);
  });
  Kokkos::fence();

  auto out_h = Kokkos::create_mirror_view(out);
  Kokkos::deep_copy(out_h, out);

  const Real tol = std::is_same<Real,double>::value ? 1e-6 : 1e-4;
  REQUIRE( std::abs(out_h(0) - static_cast<Real>(liq.hdo[i_liq])) / out_h(0) < tol );
  REQUIRE( std::abs(out_h(1) - static_cast<Real>(ice.o18[i_ice])) / out_h(1) < tol );
}

/* Verify that alternative formulations produce measurably different results.

   Deliberately pair-specific rather than an all-pairs loop over the enums: the
   bounds below are calibrated per pair, and "the formulations must differ" is
   not true in general. The two ice formulations share an identical Oxygen
   coefficient row (both take H2(18)O from Majoube 1971 and differ only in D),
   so an all-pairs "must differ" assertion would fail on O18 today. Adding a
   formulation does not require touching this function -- it just will not be
   compared here; add a block if the new pairing is worth asserting. */
void verify_formulation_differences()
{
  using Real = scream::Real;

  // Representative temperatures
  const Real t_warm = Real(293.15);  // 20°C (liquid)
  const Real t_cold = Real(243.15);  // -30°C (ice)

  /* Liquid-vapor: Horita & Wesolowski vs Majoube.

     These are two independent laboratory regressions of the same physical
     quantity, so at a temperature well inside both fitted ranges they must
     agree closely -- but not exactly, or the runtime switch would be
     meaningless. Bracketing the difference on both sides catches a silently
     ignored formulation option (lower bound) and a transcription error in
     either coefficient row (upper bound).

     The 1e-2 lower bound this assertion previously carried was calibrated when
     the Majoube path overflowed to +inf; it encoded the bug, not the physics. */
  {
    wiso::WaterIsotopeConstants<Real> const_horita;  // Default
    wiso::WaterIsotopeRuntimeOptions opts_maj;
    opts_maj.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
    wiso::WaterIsotopeConstants<Real> const_majoube(opts_maj);

    Real alpha_horita = wiso::WaterIsotopeFractionation::alpha_liquid_vapor(
      t_warm, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, const_horita);
    Real alpha_majoube = wiso::WaterIsotopeFractionation::alpha_liquid_vapor(
      t_warm, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, const_majoube);

    Real rel_diff = std::abs(alpha_horita - alpha_majoube) / alpha_horita;

    // Observed 6.2e-4 at 20 C.
    REQUIRE( rel_diff > Real(1e-5) );
    REQUIRE( rel_diff < Real(0.01) );
    // Both must still enrich the heavy isotope into the condensate.
    REQUIRE( alpha_horita > Real(1.0) );
    REQUIRE( alpha_majoube > Real(1.0) );
  }

  // Ice-vapor: Merlivat vs IsoCAM3 should differ slightly
  {
    wiso::WaterIsotopeConstants<Real> const_merlivat;  // Default
    wiso::WaterIsotopeRuntimeOptions opts_iso;
    opts_iso.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;
    wiso::WaterIsotopeConstants<Real> const_isocam3(opts_iso);

    Real alpha_merlivat = wiso::WaterIsotopeFractionation::alpha_ice_vapor(
      t_cold, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, const_merlivat);
    Real alpha_isocam3 = wiso::WaterIsotopeFractionation::alpha_ice_vapor(
      t_cold, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, const_isocam3);

    Real rel_diff = std::abs(alpha_merlivat - alpha_isocam3) / alpha_merlivat;

    // Should differ measurably but slightly (0.006% coefficient difference)
    REQUIRE( rel_diff > Real(1e-6) );
    REQUIRE( rel_diff < Real(0.01) );  // But not by more than 1%
  }
}

/* T_lane_fill/range_mask exist so that a Pack lane outside range_mask -- e.g.
   padding left by upstream physics -- can never form a 1/T term from a
   division-by-zero, which would otherwise abort an FPE-trapping build even
   though no caller reads that lane's result. Enable trapping locally (rather
   than requiring the separate SCREAM_FPE build) to prove the guard actually
   works, not just that the masked lane's numeric result looks fine. */
void verify_dead_lane_fpe_safety()
{
  using Real = scream::Real;
  using WIF  = wiso::WaterIsotopeFractionation;

  if (SCREAM_PACK_SIZE < 2) {
    return;  // no padding lane to poison
  }

  using PackN = ekat::Pack<Real, SCREAM_PACK_SIZE>;
  using MaskN = ekat::Mask<SCREAM_PACK_SIZE>;

  const int saved_fpes = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID);

  wiso::WaterIsotopeConstants<Real> constants;

  PackN t(Real(273.15));
  MaskN range_mask(true);
  t[1] = Real(0.0);          // dead lane: would force 1/T = 1/0 if not filled
  range_mask.set(1, false);

  const PackN alpha = WIF::alpha_liquid_vapor(
      t, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, constants, range_mask);

  ekat::disable_all_fpes();
  ekat::enable_fpes(saved_fpes);

  const PackN alpha_ref = WIF::alpha_liquid_vapor(
      PackN(Real(273.15)), wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor, constants);

  // Reaching this point at all (no FPE abort) is the primary check; the value
  // check confirms the poisoned dead lane did not perturb the live one.
  REQUIRE( alpha[0] == alpha_ref[0] );
}

template <typename RealT>
void run_both_pack_sizes(
  const char* phase, const AlphaTable& expected,
  RealT tol, const wiso::WaterIsotopeConstants<RealT>& constants)
{
  using Pack1 = ekat::Pack<RealT, 1>;
  using PackN = ekat::Pack<RealT, SCREAM_PACK_SIZE>;

  run_sweep<Pack1>(phase, expected, tol, constants);
  run_sweep<PackN>(phase, expected, tol, constants);
}

} // namespace

TEST_CASE("water_isotopes_fractionation") {
    using Real = scream::Real;

    /* Every sweep below compares against alpha values tabulated externally in
       fractionation_factors.xlsx, so the tolerance only needs to absorb the
       difference between two evaluation orders of the same polynomial (Horner
       here, expanded sum-of-powers in the spreadsheet), not any physics
       uncertainty. Observed worst case is ~1 ulp, so this leaves ample margin
       in both precisions. */
    const Real tight_tol = std::is_same<Real,double>::value ? Real(1e-12) : Real(1e-5);

    /* Sweep every enumerator, not a hand-written list of them, so adding a
       formulation to either enum brings it under test as soon as its expected
       values are added to the table above. */
    SECTION("all_liquid_vapor_formulations") {
      for (int f = 0; f < etoi(wiso::LiquidVaporFractionation::FormulationCount); ++f) {
        const auto formulation = static_cast<wiso::LiquidVaporFractionation>(f);
        INFO("liquid-vapor formulation index " << f << " ("
             << wiso::liquid_vapor_ref[f] << ")");
        REQUIRE( expected_is_populated(liq_expected[f]) );

        wiso::WaterIsotopeRuntimeOptions opts;
        opts.liquid_vapor = formulation;
        wiso::WaterIsotopeConstants<Real> constants(opts);
        run_both_pack_sizes("liquid-vapor", liq_expected[f], tight_tol, constants);
      }
    }

    SECTION("all_ice_vapor_formulations") {
      for (int f = 0; f < etoi(wiso::IceVaporFractionation::FormulationCount); ++f) {
        const auto formulation = static_cast<wiso::IceVaporFractionation>(f);
        INFO("ice-vapor formulation index " << f << " ("
             << wiso::ice_vapor_ref[f] << ")");
        REQUIRE( expected_is_populated(ice_expected[f]) );

        wiso::WaterIsotopeRuntimeOptions opts;
        opts.ice_vapor = formulation;
        wiso::WaterIsotopeConstants<Real> constants(opts);
        run_both_pack_sizes("ice-vapor", ice_expected[f], tight_tol, constants);
      }
    }

    /* The loops above cover the defaults only because the default enumerators
       happen to be index 0 of each enum. Assert that rather than assume it: a
       reordering that moved the default would otherwise leave the
       default-constructed path -- the one production uses -- unswept. */
    SECTION("defaults_are_covered") {
      wiso::WaterIsotopeRuntimeOptions defaults;
      REQUIRE( expected_is_populated(liq_expected[etoi(defaults.liquid_vapor)]) );
      REQUIRE( expected_is_populated(ice_expected[etoi(defaults.ice_vapor)]) );

      wiso::WaterIsotopeConstants<Real> constants;  // default-constructed
      run_both_pack_sizes("liquid-vapor",
        liq_expected[etoi(defaults.liquid_vapor)], tight_tol, constants);
      run_both_pack_sizes("ice-vapor",
        ice_expected[etoi(defaults.ice_vapor)], tight_tol, constants);
    }

    SECTION("device execution") {
      run_on_device();
    }

    SECTION("formulation_differences") {
      verify_formulation_differences();
    }

    SECTION("dead_lane_fpe_safety") {
      verify_dead_lane_fpe_safety();
    }

  }

} // namespace scream
