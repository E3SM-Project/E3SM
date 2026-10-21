#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_ocean_flux.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_pack.hpp>

#include <cmath>
#include <iomanip>
#include <iterator>
#include <limits>
#include <type_traits>
#include <vector>

namespace scream {
namespace {

using Pack1 = ekat::Pack<Real, 1>;
using WIF   = wiso::WaterIsotopeFractionation;
using wiso::WaterIsotopologues;

constexpr int num_species = etoi(WaterIsotopologues::Count);

const char* const species_names[num_species] = {"H216O", "HDO", "H218O", "H217O", "HTO"};

// Bottom-level bulk state of one ocean column, in the notation of
// wiso::isotope_ocean_flux.
struct OceanState {
  Real f, rho, u, r, s, q, t_s;
};

// Evaporating states (q < s), with ssq roughly the saturation value at the
// SST. The SSTs lie inside the fitted range of every liquid/vapor formulation,
// so the debug build prints no extrapolation warnings. q/s stays below ~0.8
// so that the cancellation in s - q does not dominate the round-off.
const OceanState evap_states[] = {
  //  f     rho   ustar  re      ssq      qv      sst [K]
  {  1.0,  1.27,  0.10, 1.0e-3, 3.8e-3,  3.0e-3, 273.15 },
  {  0.5,  1.24,  0.25, 1.2e-3, 7.6e-3,  6.0e-3, 283.15 },
  {  0.25, 1.20,  0.30, 1.5e-3, 1.45e-2, 1.0e-2, 293.15 },
  {  1.0,  1.17,  0.50, 2.0e-3, 2.2e-2,  1.7e-2, 300.15 },
  {  0.8,  1.15,  0.20, 1.8e-3, 2.9e-2,  2.2e-2, 305.0  },
};
constexpr int num_evap_states = std::size(evap_states);

// Formulations to sweep. The default gives R_ocn == 1 for every species; the
// others make R_ocn != 1, so the invariants are not checked only at R == 1.
struct Formulation {
  const char* label;
  wiso::WaterIsotopeRuntimeOptions opts;
};

std::vector<Formulation> formulations()
{
  std::vector<Formulation> out(3);
  out[0].label = "defaults (horita_wesolowski_1994, normalized, modern)";
  out[1].label = "natural_abundance + LGM";
  out[1].opts.standard_ratio   = wiso::StandardRatioFormulation::NaturalAbundance;
  out[1].opts.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
  out[2].label = "majoube_1971 + LGM";
  out[2].opts.liquid_vapor     = wiso::LiquidVaporFractionation::Majoube1971;
  out[2].opts.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
  return out;
}

// Same host-side resolution as the WaterIsotopes constructor.
Real ocean_ratio(const wiso::WaterIsotopeParameters<Real>& p, const WaterIsotopologues sp)
{
  return p.std_iso_ratio(sp) * p.mean_ocean_enrichment(sp);
}

Real bulk_evap(const OceanState& st)
{
  return water_tracers::bulk_ocean_flux(Pack1(st.f), Pack1(st.rho), Pack1(st.u),
                                        Pack1(st.r), Pack1(st.s), Pack1(st.q))[0];
}

Real alpha_lv(const OceanState& st, const WaterIsotopologues sp,
              const wiso::WaterIsotopeParameters<Real>& p)
{
  return WIF::alpha_liquid_vapor(Pack1(st.t_s), sp, wiso::CondensedOverVapor, p)[0];
}

Real iso_flux(const OceanState& st, const Real q_i, const Real e_ocn, const Real r_ocn,
              const WaterIsotopologues sp, const wiso::WaterIsotopeParameters<Real>& p)
{
  return wiso::isotope_ocean_flux(Pack1(st.f), Pack1(st.rho), Pack1(st.u), Pack1(st.r),
                                  Pack1(st.s), Pack1(st.q), Pack1(q_i), Pack1(e_ocn),
                                  Pack1(st.t_s), r_ocn, sp, p)[0];
}

// The spec's tolerance bound: |computed - expected| <= atol + rtol*|expected|.
void require_within(const Real computed, const Real expected, const Real rtol, const Real atol)
{
  const Real diff  = std::abs(computed - expected);
  const Real bound = atol + rtol*std::abs(expected);
  INFO(std::setprecision(17) << "computed = " << computed << ", expected = " << expected
       << ", |diff| = " << diff << ", bound = " << bound);
  REQUIRE(diff <= bound);
}

// Invariant (1) on device, with parameters built on the host and captured by
// value (no formulation table is indexed on device), against the host value.
void run_on_device(const wiso::WaterIsotopeRuntimeOptions& opts, const Real rtol)
{
  using KT = ekat::KokkosTypes<DefaultDevice>;

  const wiso::WaterIsotopeParameters<Real> iso_params(opts);
  Kokkos::Array<Real, num_species> r_ocn;
  for (int isp = 0; isp < num_species; ++isp) {
    r_ocn[isp] = ocean_ratio(iso_params, static_cast<WaterIsotopologues>(isp));
  }

  KT::view_1d<OceanState> states("ocean_states", num_evap_states);
  auto states_h = Kokkos::create_mirror_view(states);
  for (int i = 0; i < num_evap_states; ++i) {
    states_h(i) = evap_states[i];
  }
  Kokkos::deep_copy(states, states_h);

  KT::view_2d<Real> flux("isotope_flux", num_evap_states, num_species);
  Kokkos::parallel_for("water_isotopes_ocean_flux_device", KT::RangePolicy(0, num_evap_states),
                       KOKKOS_LAMBDA(const int i) {
    const OceanState st = states(i);
    const Pack1 e_ocn = water_tracers::bulk_ocean_flux(Pack1(st.f), Pack1(st.rho), Pack1(st.u),
                                                       Pack1(st.r), Pack1(st.s), Pack1(st.q));
    for (int isp = 0; isp < num_species; ++isp) {
      const auto sp = static_cast<WaterIsotopologues>(isp);
      const Pack1 alpha = WIF::alpha_liquid_vapor(Pack1(st.t_s), sp, wiso::CondensedOverVapor,
                                                  iso_params);
      const Pack1 q_i = (r_ocn[isp]/alpha)*st.q;
      flux(i, isp) = wiso::isotope_ocean_flux(Pack1(st.f), Pack1(st.rho), Pack1(st.u),
                                              Pack1(st.r), Pack1(st.s), Pack1(st.q), q_i, e_ocn,
                                              Pack1(st.t_s), r_ocn[isp], sp, iso_params)[0];
    }
  });
  Kokkos::fence();

  auto flux_h = Kokkos::create_mirror_view(flux);
  Kokkos::deep_copy(flux_h, flux);

  for (int i = 0; i < num_evap_states; ++i) {
    const OceanState& st = evap_states[i];
    const Real e_ocn = bulk_evap(st);
    for (int isp = 0; isp < num_species; ++isp) {
      const auto sp = static_cast<WaterIsotopologues>(isp);
      INFO("species " << species_names[isp] << ", state " << i);
      require_within(flux_h(i, isp), (r_ocn[isp]/alpha_lv(st, sp, iso_params))*e_ocn,
                     rtol, Real(0));
    }
  }
}

} // namespace

TEST_CASE("water_isotopes_ocean_flux") {
  /* The spec's bound for the equilibrium invariants is rtol = 1e-12,
     atol = 1e-15. Both sides are the same algebra evaluated in different
     orders, so only round-off separates them.

     rtol = 1e-12 is ~4500 double ulps, which single precision cannot meet, so
     an sp build gets 100 float ulps instead (about the fractionation test's
     1e-5). atol needs no such split: the invariant it applies to cancels
     exactly in either precision. */
  const Real rtol = std::is_same<Real, double>::value
                  ? Real(1e-12)
                  : Real(100)*std::numeric_limits<Real>::epsilon();
  const Real atol = Real(1e-15);

  SECTION("equilibrium_invariants") {
    for (const auto& form : formulations()) {
      INFO("formulation: " << form.label);
      const wiso::WaterIsotopeParameters<Real> iso_params(form.opts);

      for (int isp = 0; isp < num_species; ++isp) {
        const auto sp = static_cast<WaterIsotopologues>(isp);
        const Real r_ocn = ocean_ratio(iso_params, sp);
        INFO("species " << species_names[isp] << ", R_ocn = " << std::setprecision(17) << r_ocn);

        for (int i = 0; i < num_evap_states; ++i) {
          const OceanState& st = evap_states[i];
          const Real e_ocn = bulk_evap(st);
          const Real alpha = alpha_lv(st, sp, iso_params);
          INFO("state " << i << ", T_s = " << st.t_s << " K, alpha = " << alpha
               << ", E_ocn = " << e_ocn);
          REQUIRE(e_ocn > 0);

          /* (1) Vapor at the ratio of the ocean's equilibrium vapor: the flux
                 carries that ratio times the bulk evaporation. The metric is a
                 relative error, so no atol is added. At these magnitudes
                 (1e-11 to 6e-6 kg/m2/s) atol = 1e-15 would dominate the bound
                 for the rare species, so leaving it out makes the check
                 stricter. */
          const Real r_eq = r_ocn/alpha;
          require_within(iso_flux(st, r_eq*st.q, e_ocn, r_ocn, sp, iso_params),
                         r_eq*e_ocn, rtol, Real(0));

          /* (2) Vapor in equilibrium with the ocean surface: no net flux,
                 |F| <= atol. q_i is formed exactly as s_eff is inside
                 isotope_ocean_flux, so the two cancel to 0. */
          const Real q_sat_i = st.s*r_ocn/alpha;
          require_within(iso_flux(st, q_sat_i, e_ocn, r_ocn, sp, iso_params),
                         Real(0), rtol, atol);
        }
      }
    }
  }

  SECTION("h216o_matches_bulk_evaporation_exactly") {
    // Default formulations: R_ocn == 1 and alpha == 1 exactly, so the isotope
    // branch reduces to bulk_ocean_flux(f, rho, u, r, s, q) bit-for-bit.
    const wiso::WaterIsotopeParameters<Real> iso_params;
    const auto sp = WaterIsotopologues::H216O;
    const Real r_ocn = ocean_ratio(iso_params, sp);
    REQUIRE(r_ocn == Real(1));

    for (int i = 0; i < num_evap_states; ++i) {
      const OceanState& st = evap_states[i];
      INFO("state " << i);
      REQUIRE(alpha_lv(st, sp, iso_params) == Real(1));

      const Real e_ocn = bulk_evap(st);
      REQUIRE(iso_flux(st, st.q, e_ocn, r_ocn, sp, iso_params) == e_ocn);
    }
  }

  SECTION("dew_and_invalid_ocean") {
    /* Neither branch may evaluate alpha. The SST is NaN in every state below,
       so evaluating alpha would trip its temperature guard and abort. */
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    wiso::WaterIsotopeRuntimeOptions opts;
    opts.standard_ratio   = wiso::StandardRatioFormulation::NaturalAbundance;
    opts.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    const wiso::WaterIsotopeParameters<Real> iso_params(opts);

    const OceanState dew_states[] = {
      //  f    rho   ustar  re      ssq     qv      sst
      { 1.0,  1.20, 0.20, 1.5e-3, 1.0e-2, 1.1e-2, nan },  // E_ocn < 0
      { 0.3,  1.25, 0.15, 1.2e-3, 5.0e-3, 8.0e-3, nan },  // E_ocn < 0
      { 1.0,  1.20, 0.20, 1.5e-3, 1.0e-2, 1.0e-2, nan },  // E_ocn == 0
    };

    for (int isp = 0; isp < num_species; ++isp) {
      const auto sp = static_cast<WaterIsotopologues>(isp);
      const Real r_ocn = ocean_ratio(iso_params, sp);
      INFO("species " << species_names[isp]);

      for (const auto& st : dew_states) {
        const Real e_ocn = bulk_evap(st);
        INFO("dew: q = " << st.q << ", s = " << st.s << ", E_ocn = " << e_ocn);
        REQUIRE(e_ocn <= 0);

        // Dew removes the tracer at its own bottom-level ratio. Same
        // expression as the implementation, so equality is exact.
        const Real q_i = Real(0.9)*r_ocn*st.q;
        const Real expected = water_tracers::wtrc_ratio(Pack1(q_i), Pack1(st.q))[0]*e_ocn;
        REQUIRE(iso_flux(st, q_i, e_ocn, r_ocn, sp, iso_params) == expected);
      }

      // Not a valid ocean point. In the model E_ocn is also 0 here; a nonzero
      // value makes the result depend on the f gate alone.
      const OceanState land = { 0.0, 1.20, 0.20, 1.5e-3, 1.0e-2, 5.0e-3, nan };
      REQUIRE(iso_flux(land, r_ocn*land.q, Real(1e-5), r_ocn, sp, iso_params) == Real(0));
    }
  }

  SECTION("device_execution") {
    wiso::WaterIsotopeRuntimeOptions opts;
    opts.standard_ratio   = wiso::StandardRatioFormulation::NaturalAbundance;
    opts.ocean_enrichment = wiso::OceanEnrichmentFormulation::LGM;
    run_on_device(opts, rtol);
  }
}

} // namespace scream
