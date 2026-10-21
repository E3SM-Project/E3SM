#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_functions.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_fpe.hpp>
#include <ekat_pack.hpp>

#include <cfenv>
#include <limits>

namespace scream {
namespace {

namespace wt = water_tracers;
using Pack1 = ekat::Pack<Real,1>;

// An evaporating ocean column (q < s)
const Real f0 = 0.9, rho0 = 1.2, u0 = 0.3, r0 = 1.5e-3, s0 = 0.02, q0 = 0.01;

Pack1 E (const Real f, const Real q) {
  return wt::bulk_ocean_flux(Pack1(f), Pack1(rho0), Pack1(u0), Pack1(r0), Pack1(s0), Pack1(q));
}

Pack1 generic (const wt::OceanFlux kind, const Real r_ocn, const Real q, const Real q_i) {
  return wt::generic_ocean_flux(kind, r_ocn, Pack1(f0), Pack1(rho0), Pack1(u0), Pack1(r0),
                                Pack1(s0), Pack1(q), Pack1(q_i), E(f0, q));
}

TEST_CASE("generic_ocean_flux") {
  using OF = wt::OceanFlux;

  SECTION("bulk_formula with R_ocn=1, q_i=q equals E_ocn bit-for-bit") {
    REQUIRE(generic(OF::BulkFormula, 1, q0, q0)[0] == E(f0, q0)[0]);
  }
  SECTION("bulk_proportional is R_ocn*E_ocn") {
    REQUIRE(generic(OF::BulkProportional, 0.5, q0, 3e-3)[0] == Real(0.5)*E(f0, q0)[0]);
    REQUIRE(generic(OF::BulkProportional, 1,   q0, 3e-3)[0] == E(f0, q0)[0]);
  }
  SECTION("bulk_formula uses s*R_ocn and q_i") {
    const Real r_ocn = 0.5, q_i = 2e-3;
    const Real ref = wt::bulk_ocean_flux(Pack1(f0), Pack1(rho0), Pack1(u0), Pack1(r0),
                                         Pack1(s0*r_ocn), Pack1(q_i))[0];
    REQUIRE(generic(OF::BulkFormula, r_ocn, q0, q_i)[0] == ref);
  }
  SECTION("dew removes the tracer at its ratio to qv") {
    const Real q = 0.03, q_i = 0.012;  // q > s: dew
    const Real e = E(f0, q)[0];
    REQUIRE(e < 0);
    const Real ref = (wt::wtrc_ratio(Pack1(q_i), Pack1(q))*e)[0];
    REQUIRE(generic(OF::BulkProportional, 0.5, q, q_i)[0] == ref);
    REQUIRE(generic(OF::BulkFormula,      0.5, q, q_i)[0] == ref);
  }
}

TEST_CASE("invalid_ocean") {
  // Trap explicitly: no step may compute with fill values or NaN.
  const int saved = ekat::get_enabled_fpes();
  ekat::disable_all_fpes();
  ekat::enable_fpes(FE_DIVBYZERO | FE_INVALID);

  const Real nan  = std::numeric_limits<Real>::quiet_NaN();
  const Real fill = 1e30;
  struct State { Real f, u, r, s; };
  const State states[] = {
    {0,  u0, r0,   s0},    // no ocean
    {f0, 0,  r0,   s0},    // no friction velocity
    {f0, u0, r0,   fill},  // coupler fill value
    {f0, u0, fill, s0},
    {f0, u0, nan,  s0},    // NaN
    {f0, u0, r0,   nan},
    {nan, u0, r0,  s0},
  };

  for (const auto& st : states) {
    Pack1 f(st.f), u(st.u), r(st.r), s(st.s);
    REQUIRE(not wt::ocean_valid(f, u, r, s)[0]);

    const Pack1 e = wt::ocean_evaporation(f, u, r, s, Pack1(rho0), Pack1(q0));
    REQUIRE(e[0] == 0);
    REQUIRE(f[0] == 0);
    REQUIRE(u[0] == 0);
    REQUIRE(r[0] == 0);
    REQUIRE(s[0] == 0);

    for (const auto kind : {wt::OceanFlux::BulkProportional, wt::OceanFlux::BulkFormula}) {
      REQUIRE(wt::generic_ocean_flux(kind, Real(1), f, Pack1(rho0), u, r, s, Pack1(q0), Pack1(q0), e)[0] == 0);
    }
  }

  // A valid state is left untouched, and E_ocn is the bulk formula exactly
  Pack1 f(f0), u(u0), r(r0), s(s0);
  REQUIRE(wt::ocean_valid(f, u, r, s)[0]);
  const Pack1 e = wt::ocean_evaporation(f, u, r, s, Pack1(rho0), Pack1(q0));
  REQUIRE(f[0] == f0);
  REQUIRE(e[0] == E(f0, q0)[0]);

  ekat::disable_all_fpes();
  ekat::enable_fpes(saved);
}

} // anonymous namespace
} // namespace scream
