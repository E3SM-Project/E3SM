#include "catch2/catch.hpp"

#include "physics/specialized_tracers/water_tracers/eamxx_water_tracers_functions.hpp"
#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_fractionation.hpp"
#include "physics/specialized_tracers/water_isotopes/eamxx_water_isotopes_parameters.hpp"

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
using WIF   = water_isotopes::WaterIsotopeFractionation;
using water_isotopes::WaterIsotopologues;

constexpr int num_species = etoi(WaterIsotopologues::Count);

const char* const species_names[num_species] = {"H216O", "HDO", "H218O", "H217O", "HTO"};

// Bottom-level bulk state of one ocean column, in the notation of
// water_tracers::ocean_flux.
struct OceanState {
  Real f, rho, u, r, s, q, t_s;
};

// Evaporating states (q < s), with ssq roughly the saturation value at the
// SST. The first is at the sea-water freezing point, below the fitted range
// of every liquid/vapor formulation but above its Tmin_warn, so the debug
// build prints no extrapolation warnings. q/s stays below ~0.8 so that the
// cancellation in s - q does not dominate the round-off.
const OceanState evap_states[] = {
  //  f     rho   ustar  re      ssq      qv      sst [K]
  {  1.0,  1.29,  0.10, 1.0e-3, 3.2e-3,  2.5e-3, 271.35 },
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
  water_isotopes::WaterIsotopeRuntimeOptions opts;
};

std::vector<Formulation> formulations()
{
  std::vector<Formulation> out(3);
  out[0].label = "defaults (horita_wesolowski_1994, normalized, modern)";
  out[1].label = "natural_abundance + LGM";
  out[1].opts.standard_ratio   = water_isotopes::StandardRatioFormulation::NaturalAbundance;
  out[1].opts.ocean_enrichment = water_isotopes::OceanEnrichmentFormulation::LGM;
  out[2].label = "majoube_1971 + LGM";
  out[2].opts.liquid_vapor     = water_isotopes::LiquidVaporFractionation::Majoube1971;
  out[2].opts.ocean_enrichment = water_isotopes::OceanEnrichmentFormulation::LGM;
  return out;
}

// Same host-side resolution as the WaterIsotopes constructor.
Real ocean_ratio(const water_isotopes::WaterIsotopeParameters<Real>& p, const WaterIsotopologues sp)
{
  return p.std_iso_ratio(sp) * p.mean_ocean_enrichment(sp);
}

Real bulk_evap(const OceanState& st)
{
  return water_tracers::bulk_ocean_flux(Pack1(st.f), Pack1(st.rho), Pack1(st.u),
                                        Pack1(st.r), Pack1(st.s), Pack1(st.q))[0];
}

Real alpha_lv(const OceanState& st, const WaterIsotopologues sp,
              const water_isotopes::WaterIsotopeParameters<Real>& p)
{
  return WIF::alpha_liquid_vapor(Pack1(st.t_s), sp, water_isotopes::CondensedOverVapor, p)[0];
}

Real iso_flux(const OceanState& st, const Real q_i, const Real e_ocn, const Real r_ocn,
              const WaterIsotopologues sp, const water_isotopes::WaterIsotopeParameters<Real>& p)
{
  const Pack1 f(st.f), e(e_ocn);
  const auto evap = water_tracers::is_evaporating(f, e);
  const Pack1 alpha = WIF::alpha_liquid_vapor(Pack1(st.t_s), sp, water_isotopes::CondensedOverVapor,
                                              p, evap);
  return water_tracers::ocean_flux(water_tracers::OceanFlux::BulkFormula, r_ocn,
      f, Pack1(st.rho), Pack1(st.u), Pack1(st.r), Pack1(st.s), Pack1(st.q),
      Pack1(q_i), e, alpha)[0];
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
void run_on_device(const water_isotopes::WaterIsotopeRuntimeOptions& opts, const Real rtol)
{
  using KT = ekat::KokkosTypes<DefaultDevice>;

  const water_isotopes::WaterIsotopeParameters<Real> iso_params(opts);
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
    const Pack1 f(st.f);
    const Pack1 e_ocn = water_tracers::bulk_ocean_flux(f, Pack1(st.rho), Pack1(st.u),
                                                       Pack1(st.r), Pack1(st.s), Pack1(st.q));
    const auto evap = water_tracers::is_evaporating(f, e_ocn);
    for (int isp = 0; isp < num_species; ++isp) {
      const auto sp = static_cast<WaterIsotopologues>(isp);
      const Pack1 alpha = WIF::alpha_liquid_vapor(Pack1(st.t_s), sp, water_isotopes::CondensedOverVapor,
                                                  iso_params, evap);
      const Pack1 q_i = (r_ocn[isp]/alpha)*st.q;
      flux(i, isp) = water_tracers::ocean_flux(water_tracers::OceanFlux::BulkFormula, r_ocn[isp],
          f, Pack1(st.rho), Pack1(st.u), Pack1(st.r), Pack1(st.s), Pack1(st.q), q_i, e_ocn,
          alpha)[0];
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
      const water_isotopes::WaterIsotopeParameters<Real> iso_params(form.opts);

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
                 ocean_flux, so the two cancel to 0. */
          const Real q_sat_i = st.s*r_ocn/alpha;
          require_within(iso_flux(st, q_sat_i, e_ocn, r_ocn, sp, iso_params),
                         Real(0), rtol, atol);
        }
      }
    }
  }

  SECTION("device_execution") {
    water_isotopes::WaterIsotopeRuntimeOptions opts;
    opts.standard_ratio   = water_isotopes::StandardRatioFormulation::NaturalAbundance;
    opts.ocean_enrichment = water_isotopes::OceanEnrichmentFormulation::LGM;
    run_on_device(opts, rtol);
  }
}

} // namespace scream
