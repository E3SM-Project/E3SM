#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_fractionation.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_pack.hpp>
#include <ekat_view_utils.hpp>

#include <cmath>

namespace scream {
namespace {

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

// Temperature sweeps.
constexpr int NLIQ = 6;
constexpr int NICE = 6;
const double T_liq[NLIQ] = {233.15, 253.15, 273.15, 283.15, 293.15, 303.15};
const double T_ice[NICE] = {213.15, 233.15, 243.15, 253.15, 263.15, 273.15};

template <typename ScalarT>
void run_sweep()
{
  using WIF     = wiso::WaterIsotopeFractionation;
  using STraits = ekat::ScalarTraits<ScalarT>;
  using RealT   = typename STraits::scalar_type;

  // Relative tolerance: tighter in double precision, looser in single.
  const RealT tol = std::is_same<RealT,double>::value ? 1e-6 : 1e-4;

  auto rel_ok = [&](const ScalarT& computed, double expected) {
    // ScalarT is a Pack of identical lanes (all set from one temperature),
    // so comparing lane 0 is sufficient and works for Real too.
    const RealT c = computed[0];
    return std::abs(c - static_cast<RealT>(expected)) / std::abs(static_cast<RealT>(expected)) < tol;
  };

  using wiso::CondensedOverVapor;
  using wiso::VaporOverCondensed;
  using wiso::H216O;
  using wiso::HDO;
  using wiso::H218O;
  using wiso::H217O;
  using wiso::HTO;

  // ---- Liquid-vapor sweep ----
  double prev_hdo = 1e30, prev_o18 = 1e30;
  for (int i = 0; i < NLIQ; ++i) {
    const ScalarT t(T_liq[i]);

    const ScalarT a_hdo = WIF::alpha_liquid_vapor(t, HDO,   CondensedOverVapor);
    const ScalarT a_o18 = WIF::alpha_liquid_vapor(t, H218O, CondensedOverVapor);

    // Absolute-value checks against the independent reference.
    REQUIRE( rel_ok(a_hdo, ref_alpl_hdo(T_liq[i])) );
    REQUIRE( rel_ok(a_o18, ref_alpl_o18(T_liq[i])) );

    // Ordinary water is non-fractionating, exactly.
    const ScalarT a_16 = WIF::alpha_liquid_vapor(t, H216O, CondensedOverVapor);
    REQUIRE( (a_16 == ScalarT(1)).all() );

    // Direction: VaporOverCondensed is the reciprocal.
    const ScalarT a_hdo_inv = WIF::alpha_liquid_vapor(t, HDO,   VaporOverCondensed);
    const ScalarT a_o18_inv = WIF::alpha_liquid_vapor(t, H218O, VaporOverCondensed);
    REQUIRE( rel_ok(a_hdo_inv, 1.0/ref_alpl_hdo(T_liq[i])) );
    REQUIRE( rel_ok(a_o18_inv, 1.0/ref_alpl_o18(T_liq[i])) );

    // Derived species power laws.
    const ScalarT a_17 = WIF::alpha_liquid_vapor(t, H217O, CondensedOverVapor);
    const ScalarT a_ht = WIF::alpha_liquid_vapor(t, HTO,   CondensedOverVapor);
    REQUIRE( rel_ok(a_17, std::pow(ref_alpl_o18(T_liq[i]), 0.529)) );
    REQUIRE( rel_ok(a_ht, std::pow(ref_alpl_hdo(T_liq[i]), 2.0)) );

    // Monotonic decrease toward 1 with increasing T, and alpha >= 1.
    REQUIRE( a_hdo[0] >= RealT(1) );
    REQUIRE( a_o18[0] >= RealT(1) );
    REQUIRE( a_hdo[0] < prev_hdo );
    REQUIRE( a_o18[0] < prev_o18 );
    prev_hdo = a_hdo[0];
    prev_o18 = a_o18[0];
  }

  // ---- Ice-vapor sweep ----
  prev_hdo = 1e30; prev_o18 = 1e30;
  for (int i = 0; i < NICE; ++i) {
    const ScalarT t(T_ice[i]);

    const ScalarT a_hdo = WIF::alpha_ice_vapor(t, HDO,   CondensedOverVapor);
    const ScalarT a_o18 = WIF::alpha_ice_vapor(t, H218O, CondensedOverVapor);

    REQUIRE( rel_ok(a_hdo, ref_alpi_hdo(T_ice[i])) );
    REQUIRE( rel_ok(a_o18, ref_alpi_o18(T_ice[i])) );

    const ScalarT a_16 = WIF::alpha_ice_vapor(t, H216O, CondensedOverVapor);
    REQUIRE( (a_16 == ScalarT(1)).all() );

    const ScalarT a_hdo_inv = WIF::alpha_ice_vapor(t, HDO,   VaporOverCondensed);
    const ScalarT a_o18_inv = WIF::alpha_ice_vapor(t, H218O, VaporOverCondensed);
    REQUIRE( rel_ok(a_hdo_inv, 1.0/ref_alpi_hdo(T_ice[i])) );
    REQUIRE( rel_ok(a_o18_inv, 1.0/ref_alpi_o18(T_ice[i])) );

    const ScalarT a_17 = WIF::alpha_ice_vapor(t, H217O, CondensedOverVapor);
    const ScalarT a_ht = WIF::alpha_ice_vapor(t, HTO,   CondensedOverVapor);
    REQUIRE( rel_ok(a_17, std::pow(ref_alpi_o18(T_ice[i]), 0.529)) );
    REQUIRE( rel_ok(a_ht, std::pow(ref_alpi_hdo(T_ice[i]), 2.0)) );

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
  Kokkos::parallel_for("wiso_frac_device", 1, KOKKOS_LAMBDA(const int /*i*/) {
    out(0) = WIF::alpha_liquid_vapor(Real(273.15), wiso::HDO,   wiso::CondensedOverVapor);
    out(1) = WIF::alpha_ice_vapor  (Real(253.15), wiso::H218O, wiso::CondensedOverVapor);
  });
  Kokkos::fence();

  auto out_h = Kokkos::create_mirror_view(out);
  Kokkos::deep_copy(out_h, out);

  const Real tol = std::is_same<Real,double>::value ? 1e-6 : 1e-4;
  REQUIRE( std::abs(out_h(0) - static_cast<Real>(ref_alpl_hdo(273.15))) / out_h(0) < tol );
  REQUIRE( std::abs(out_h(1) - static_cast<Real>(ref_alpi_o18(253.15))) / out_h(1) < tol );
}

} // anonymous namespace

TEST_CASE("water_isotopes_fractionation", "[water_isotopes]")
{
  // Real (as a length-1 pack view via operator[]) and Pack instantiations.
  run_sweep<ekat::Pack<Real,1>>();
  run_sweep<ekat::Pack<Real,SCREAM_PACK_SIZE>>();
  run_on_device();
}

} // namespace scream
