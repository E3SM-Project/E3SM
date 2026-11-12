#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_functions.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_fpe.hpp>
#include <ekat_kokkos_types.hpp>
#include <ekat_pack.hpp>

#include <cfenv>
#include <limits>

namespace scream {
namespace {

namespace wt = water_tracers;
using Pack1 = ekat::Pack<Real,1>;
using Pack4 = ekat::Pack<Real,4>;

// Enables FE_DIVBYZERO|FE_INVALID for its lifetime (not just in the fpe
// build), so the masked-division tests prove they raise nothing.
struct FpeTrap {
  int saved;
  FpeTrap () : saved(ekat::get_enabled_fpes()) {
    ekat::disable_all_fpes();
    ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID);
  }
  ~FpeTrap () {
    ekat::disable_all_fpes();
    ekat::enable_fpes(saved);
  }
};

TEST_CASE("tracer_ratio") {
  SECTION("normal ratio") {
    const Pack1 r = wt::tracer_ratio(Pack1(2e-3), Pack1(4e-3));
    REQUIRE(r[0] == Real(2e-3)/Real(4e-3));
  }
  SECTION("tiny denominator returns r_std") {
    REQUIRE(wt::tracer_ratio(Pack1(1e-24), Pack1(1e-23))[0] == Real(1));
    REQUIRE(wt::tracer_ratio(Pack1(1e-24), Pack1(-1e-23))[0] == Real(1));
    REQUIRE(wt::tracer_ratio(Pack1(1e-24), Pack1(1e-23), Real(0.5))[0] == Real(0.5));
  }
  SECTION("dry cell is FPE-safe") {
    FpeTrap trap;
    REQUIRE(wt::tracer_ratio(Pack1(0), Pack1(0))[0] == Real(1));
    REQUIRE(wt::tracer_ratio(Pack1(0), Pack1(0), Real(0.25))[0] == Real(0.25));
  }
  SECTION("mixed lanes") {
    FpeTrap trap;
    Pack4 qtrc, qtot;
    qtrc[0] = 1e-3; qtot[0] = 2e-3;   // normal
    qtrc[1] = 0;    qtot[1] = 0;      // dry: 0/0 if not masked
    qtrc[2] = 5e-4; qtot[2] = 1e-23;  // below water_tracers_qmin
    qtrc[3] = 3e-3; qtot[3] = 3e-3;   // ratio 1
    const Pack4 r = wt::tracer_ratio(qtrc, qtot, Real(0.7));
    REQUIRE(r[0] == Real(1e-3)/Real(2e-3));
    REQUIRE(r[1] == Real(0.7));
    REQUIRE(r[2] == Real(0.7));
    REQUIRE(r[3] == Real(1));
  }
}

TEST_CASE("bulk_ocean_flux") {
  const Real f = 0.8, rho = 1.2, u = 0.3, r = 1.5e-3, s = 0.02, q = 0.012;
  const Pack1 e = wt::bulk_ocean_flux(Pack1(f), Pack1(rho), Pack1(u), Pack1(r), Pack1(s), Pack1(q));

  // Upward positive, in the documented operation order
  REQUIRE(e[0] == f*rho*u*r*(s - q));
  REQUIRE(e[0] > 0);

  // s_eff = s*1, q_i = q gives E_ocn bit-for-bit
  const Pack1 e1 = wt::bulk_ocean_flux(Pack1(f), Pack1(rho), Pack1(u), Pack1(r), Pack1(s*Real(1)), Pack1(q));
  REQUIRE(e1[0] == e[0]);

  // Equilibrium: no flux when q_i = s_eff
  REQUIRE(wt::bulk_ocean_flux(Pack1(f), Pack1(rho), Pack1(u), Pack1(r), Pack1(s), Pack1(s))[0] == 0);

  // Dew when the air is above saturation
  REQUIRE(wt::bulk_ocean_flux(Pack1(f), Pack1(rho), Pack1(u), Pack1(r), Pack1(q), Pack1(s))[0] < 0);
}

TEST_CASE("functions_on_device") {
  using KT = ekat::KokkosTypes<DefaultDevice>;

  // Lane 0: normal; lane 1: dry cell; lane 2: E_ocn identity
  KT::view_1d<Real> out("out", 3);
  Kokkos::parallel_for(KT::RangePolicy(0, 1), KOKKOS_LAMBDA(const int) {
    out(0) = wt::tracer_ratio(Pack1(1e-3), Pack1(4e-3))[0];
    out(1) = wt::tracer_ratio(Pack1(0), Pack1(0), Real(0.5))[0];
    const Pack1 e  = wt::bulk_ocean_flux(Pack1(1), Pack1(1.2), Pack1(0.3), Pack1(1.5e-3), Pack1(0.02), Pack1(0.01));
    const Pack1 e1 = wt::bulk_ocean_flux(Pack1(1), Pack1(1.2), Pack1(0.3), Pack1(1.5e-3), Pack1(Real(0.02)*Real(1)), Pack1(0.01));
    out(2) = (e1[0] == e[0]) ? 1 : 0;
  });
  const auto out_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), out);
  REQUIRE(out_h(0) == Real(1e-3)/Real(4e-3));
  REQUIRE(out_h(1) == Real(0.5));
  REQUIRE(out_h(2) == 1);
}

} // anonymous namespace
} // namespace scream
