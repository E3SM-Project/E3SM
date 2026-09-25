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
  using RealT = typename ekat::ScalarTraits<ScalarT>::scalar_type;
  const auto expected_pack = ScalarT(static_cast<RealT>(expected));
  const auto rel_err = ekat::abs((computed - expected_pack) / expected_pack);
  return (rel_err < tol).all();  // Check ALL lanes, not just [0]
} 

// Independent reference implementations (formulas transcribed independently
// from the ported code) so this test doubles as a coefficient-transcription
// guard. T in Kelvin; return the raw table value R_condensed/R_vapor (>= 1).
double ref_alpl_hdo(double t) {
  return std::exp(1158.8e-12*t*t*t - 1620.1e-9*t*t + 794.84e-6*t - 161.04e-3
                  + 2.9992e6/(t*t*t));
}
double ref_alpl_o18(double t) {
  return std::exp(0.35041e6/(t*t*t) - 1.6664e3/(t*t) + 6.7123/t - 7.685e-3);
}
double ref_alpi_hdo(double t) {
  return std::exp(16289.0/(t*t) - 9.45e-2);
}
double ref_alpi_o18(double t) {
  return std::exp(11.839/t - 28.224e-3);
}

// Majoube 1971 reference implementations (alternative liquid-vapor formulation).
// alpha = exp(A/T^2 + B/T + C)
double ref_alpl_hdo_majoube(double t) {
  return std::exp(24.844e3/(t*t) - 76.248/t + 52.612e-3);
}
double ref_alpl_o18_majoube(double t) {
  return std::exp(1.137e3/(t*t) - 0.4156/t - 2.0667e-3);
}

// IsoCAM3 reference implementations (alternative ice-vapor formulation).
// alpha = exp(A/T^2 + B/T + C)
double ref_alpi_hdo_isocam3(double t) {
  return std::exp(16288.0/(t*t) - 9.34e-2);
}
double ref_alpi_o18_isocam3(double t) {
  return std::exp(11.839/t - 28.224e-3);  // Same as Merlivat for O18
}

// Temperature sweeps.
constexpr int NLIQ = 6;
constexpr int NICE = 6;
const double T_liq[NLIQ] = {233.15, 253.15, 273.15, 283.15, 293.15, 303.15};
const double T_ice[NICE] = {213.15, 233.15, 243.15, 253.15, 263.15, 273.15};

template <typename ScalarT>
void run_sweep(
  const char* phase_name,  // "liquid-vapor" or "ice-vapor"
  const double* T_array,    // Temperature array
  int N,                    // Array length
  std::function<double(double)> ref_hdo,   // HDO reference function
  std::function<double(double)> ref_o18,   // O18 reference function
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

  double prev_hdo = 1e30, prev_o18 = 1e30;
  for (int i = 0; i < N; ++i) {
    const ScalarT t(T_array[i]);

    const ScalarT a_hdo = alpha_fn(t, wiso::WaterIsotopologues::HDO, wiso::CondensedOverVapor,
constants);
    const ScalarT a_o18 = alpha_fn(t, wiso::WaterIsotopologues::H218O, wiso::CondensedOverVapor,
constants);

    REQUIRE( relative_approx(a_hdo, ref_hdo(T_array[i]), tol) );
    REQUIRE( relative_approx(a_o18, ref_o18(T_array[i]), tol) );
    
    // H216O check
    const ScalarT a_16 = alpha_fn(t, wiso::WaterIsotopologues::H216O, wiso::CondensedOverVapor,
constants);
    REQUIRE( (a_16 == ScalarT(1)).all() );
    
    // Direction checks
    const ScalarT a_hdo_inv = alpha_fn(t, wiso::WaterIsotopologues::HDO, wiso::VaporOverCondensed,
constants);
    const ScalarT a_o18_inv = alpha_fn(t, wiso::WaterIsotopologues::H218O, wiso::VaporOverCondensed,
constants);
    REQUIRE( relative_approx(a_hdo_inv, 1.0/ref_hdo(T_array[i]), tol) );
    REQUIRE( relative_approx(a_o18_inv, 1.0/ref_o18(T_array[i]), tol) );
    
    // Power law checks for H217O and HTO
    const ScalarT a_17 = alpha_fn(t, wiso::WaterIsotopologues::H217O, wiso::CondensedOverVapor,
constants);
    const ScalarT a_ht = alpha_fn(t, wiso::WaterIsotopologues::HTO, wiso::CondensedOverVapor,
constants);
    REQUIRE( relative_approx(a_17, std::pow(ref_o18(T_array[i]), 0.529), tol) );
    REQUIRE( relative_approx(a_ht, std::pow(ref_hdo(T_array[i]), 2.0), tol) );

    // Monotonicity and >= 1 checks
    REQUIRE( a_hdo[0] >= RealT(1) );
    REQUIRE( a_o18[0] >= RealT(1) );
    REQUIRE( a_hdo[0] < prev_hdo );
    REQUIRE( a_o18[0] < prev_o18 );
    prev_hdo = a_hdo[0];
    prev_o18 = a_o18[0];
  } 
} 
// Exercise the KOKKOS_INLINE_FUNCTION on the device to confirm it is
// device-callable and gives the same result as the host path.
void run_on_device()
{
  using WIF = wiso::WaterIsotopeFractionation;
  using KT  = ekat::KokkosTypes<DefaultDevice>;
  using view_1d = typename KT::template view_1d<Real>;

  view_1d out("wiso_alpha_device", 2);
  // Constructed inside the kernel to confirm the resolving constructor itself is
  // device-callable, not just the evaluator.
  Kokkos::parallel_for("wiso_frac_device", 1, KOKKOS_LAMBDA(const int /*i*/) {
    wiso::WaterIsotopeConstants<Real> constants;
    out(0) = WIF::alpha_liquid_vapor(Real(273.15), wiso::WaterIsotopologues::HDO,   wiso::CondensedOverVapor, constants);
    out(1) = WIF::alpha_ice_vapor  (Real(253.15), wiso::WaterIsotopologues::H218O, wiso::CondensedOverVapor, constants);
  });
  Kokkos::fence();

  auto out_h = Kokkos::create_mirror_view(out);
  Kokkos::deep_copy(out_h, out);

  const Real tol = std::is_same<Real,double>::value ? 1e-6 : 1e-4;
  REQUIRE( std::abs(out_h(0) - static_cast<Real>(ref_alpl_hdo(273.15))) / out_h(0) < tol );
  REQUIRE( std::abs(out_h(1) - static_cast<Real>(ref_alpi_o18(253.15))) / out_h(1) < tol );
}

// Verify that alternative formulations produce measurably different results.
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
  const char* phase, const double* T_array, int N,
  std::function<double(double)> ref_hdo, std::function<double(double)> ref_o18,
  RealT tol, const wiso::WaterIsotopeConstants<RealT>& constants)
{ 
  using Pack1 = ekat::Pack<RealT, 1>;
  using PackN = ekat::Pack<RealT, SCREAM_PACK_SIZE>;
  
  run_sweep<Pack1>(phase, T_array, N, ref_hdo, ref_o18, tol, constants);
  run_sweep<PackN>(phase, T_array, N, ref_hdo, ref_o18, tol, constants);
}

} // namespace

TEST_CASE("water_isotopes_fractionation") {
    using Real = scream::Real;

    /* Every sweep below compares against an independently transcribed reference,
       so the tolerance only needs to absorb the difference between two
       evaluation orders of the same polynomial (Horner here, expanded in the
       reference), not any physics uncertainty. Observed worst case is ~1.4 ulp,
       so this leaves ample margin in both precisions. */
    const Real tight_tol = std::is_same<Real,double>::value ? Real(1e-12) : Real(1e-5);

    SECTION("default_formulations") {
      wiso::WaterIsotopeConstants<Real> constants;
      run_both_pack_sizes("liquid-vapor", T_liq, NLIQ, ref_alpl_hdo, ref_alpl_o18,
        tight_tol, constants);
      run_both_pack_sizes("ice-vapor", T_ice, NICE, ref_alpi_hdo, ref_alpi_o18,
        tight_tol, constants);
    }

    SECTION("alternative_liquid_vapor") {
      wiso::WaterIsotopeRuntimeOptions opts;
      opts.liquid_vapor = wiso::LiquidVaporFractionation::Majoube1971;
      wiso::WaterIsotopeConstants<Real> constants(opts);
      run_both_pack_sizes("liquid-vapor", T_liq, NLIQ, ref_alpl_hdo_majoube,
        ref_alpl_o18_majoube, tight_tol, constants);
    }
    
    SECTION("alternative_ice_vapor") {
      wiso::WaterIsotopeRuntimeOptions opts;
      opts.ice_vapor = wiso::IceVaporFractionation::IsoCAM3;
      wiso::WaterIsotopeConstants<Real> constants(opts);
      run_both_pack_sizes("ice-vapor", T_ice, NICE, ref_alpi_hdo_isocam3,
        ref_alpi_o18_isocam3, tight_tol, constants);
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
