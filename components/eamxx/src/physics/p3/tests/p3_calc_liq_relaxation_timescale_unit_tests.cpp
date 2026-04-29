#include "catch2/catch.hpp"

#include "p3_functions.hpp"
#include "p3_test_data.hpp"
#include "p3_unit_tests_common.hpp"

#include "share/core/eamxx_types.hpp"

#include <cmath>
#include <limits>
#include <array>
#include <algorithm>

/*
 * Property tests for calc_liq_relaxation_timescale.
 *
 * For active rain lanes,
 *
 *   epsr = 2*pi*cdistr*rho*dv *
 *          ( f1r * Gamma(mu_r + 2)/lamr
 *          + f2r * sqrt(rho/mu) * cbrt(sc) * T(mu_r, lamr) )
 *
 * where T(mu_r, lamr) is the rain evaporation lookup-table value.
 *
 * For active cloud-liquid lanes,
 *
 *   epsc = 2*pi*rho*dv*cdist.
 *
 * These tests check algebraic identities and physical bounds implied by
 * those formulas. They intentionally keep the production kernel as the
 * source of truth and compare against independent identities, scalings,
 * and threshold behavior.
 */

namespace scream {
namespace p3 {
namespace unit_test {

namespace {

constexpr Real kAbsoluteTol = 1e-30;

template <typename T>
KOKKOS_INLINE_FUNCTION T max3(const T& a, const T& b, const T& c) {
  return a > b ? (a > c ? a : c) : (b > c ? b : c);
}

} // anonymous namespace

template <typename D>
struct UnitWrap::UnitTest<D>::TestCalcLiqRelaxationTimescale : public UnitWrap::UnitTest<D>::Base {

  struct RelaxationLane {
    Real rho, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld;
    bool context;
    Real epsr_in, epsc_in;
    Real epsr_out, epsc_out;
  };

  // Build one deterministic pack of physically valid active lanes.
  //
  // Each lane uses slightly different positive values. This avoids a test
  // suite that only validates one scalar point, while still keeping the
  // expected algebra simple. qr_incld and qc_incld start above QSMALL, and
  // context starts true, so each lane is active unless a SECTION modifies it.
  std::array<RelaxationLane, Pack::n> make_active_lanes() const
  {
    std::array<RelaxationLane, Pack::n> lanes;
    for (Int s = 0; s < Pack::n; ++s) {
      lanes[s] = {
        0.95 + 0.07 * s,
        1.2e-5 * (1 + 0.1 * s),
        1.6e-5 * (1 + 0.08 * s),
        0.75 + 0.05 * s,
        1.4 + 0.1 * s,
        800.0 + 25.0 * s,
        0.85 + 0.04 * s,
        0.65 + 0.03 * s,
        5.0 * C::QSMALL,
        4.0 * C::QSMALL,
        true,
        -999.0,
        -12345.0,
        0.0,
        0.0
      };
    }
    return lanes;
  }

  void run_kernel(const view_2d_table& revap_table_vals, const Scalar& f1r,
                  const Scalar& f2r, std::array<RelaxationLane, Pack::n>& lanes) const
  {
    // Run the production kernel on one pack and copy epsr/epsc back into lanes.
    //
    // The property tests modify lane inputs on the host, call this helper, and
    // then check the outputs. Running through Kokkos keeps the test close to the
    // real packed implementation instead of reimplementing the routine on host.
    using KTH = KokkosTypes<HostDevice>;

    KTH::view_1d<RelaxationLane> lanes_host("lanes_host", Pack::n);
    view_1d<RelaxationLane> lanes_device("lanes_device", Pack::n);
    std::copy(lanes.begin(), lanes.end(), lanes_host.data());
    Kokkos::deep_copy(lanes_device, lanes_host);

    Kokkos::parallel_for(1, KOKKOS_LAMBDA(const Int&) {
      Pack rho, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld;
      Pack epsr, epsc;
      Mask context;

      for (Int s = 0; s < Pack::n; ++s) {
        rho[s]      = lanes_device(s).rho;
        dv[s]       = lanes_device(s).dv;
        mu[s]       = lanes_device(s).mu;
        sc[s]       = lanes_device(s).sc;
        mu_r[s]     = lanes_device(s).mu_r;
        lamr[s]     = lanes_device(s).lamr;
        cdistr[s]   = lanes_device(s).cdistr;
        cdist[s]    = lanes_device(s).cdist;
        qr_incld[s] = lanes_device(s).qr_incld;
        qc_incld[s] = lanes_device(s).qc_incld;
        epsr[s]     = lanes_device(s).epsr_in;
        epsc[s]     = lanes_device(s).epsc_in;
        context.set(s, lanes_device(s).context);
      }

      Functions::calc_liq_relaxation_timescale(revap_table_vals, rho, f1r, f2r, dv,
                                               mu, sc, mu_r, lamr, cdistr, cdist,
                                               qr_incld, qc_incld, epsr, epsc,
                                               context);

      for (Int s = 0; s < Pack::n; ++s) {
        lanes_device(s).epsr_out = epsr[s];
        lanes_device(s).epsc_out = epsc[s];
      }
    });
    Kokkos::deep_copy(lanes_host, lanes_device);
    std::copy(lanes_host.data(), lanes_host.data() + Pack::n, lanes.begin());
  }

  void run_phys()
  {
    const auto& revap_table_vals = this->lookup_tables.revap_table_vals;

    // Algebraic identities should hold to roundoff. kAbsoluteTol protects
    // exact-zero checks in near-zero regimes, while identity_tol scales with
    // the magnitude of the compared values.
    const Real identity_tol = 100 * std::numeric_limits<Real>::epsilon();

    // Inequality checks need a small tolerance so that roundoff does not
    // turn an exact bound into a false failure.
    const Real inequality_tol = 100 * std::numeric_limits<Real>::epsilon();

    const auto near_tol = [&](const Real actual, const Real expected) {
      return kAbsoluteTol + identity_tol * max3<Real>(1, std::abs(actual), std::abs(expected));
    };

    const auto require_near = [&](const Real actual, const Real expected) {
      INFO("actual=" << actual << " expected=" << expected
           << " tol=" << near_tol(actual, expected));
      REQUIRE(std::abs(actual - expected) <= near_tol(actual, expected));
    };

    const auto require_zero = [&](const Real actual) {
      INFO("actual=" << actual << " abs_tol=" << kAbsoluteTol);
      REQUIRE(std::abs(actual) <= kAbsoluteTol);
    };

    const auto require_nonnegative = [&](const Real actual) {
      INFO("actual=" << actual << " tol=" << inequality_tol);
      REQUIRE(actual >= -inequality_tol);
    };

    const auto require_finite = [&](const Real actual) {
      INFO("actual=" << actual);
      REQUIRE(std::isfinite(actual));
    };

    const auto require_le = [&](const Real lhs, const Real rhs) {
      const Real tol = kAbsoluteTol + inequality_tol * max3<Real>(1, std::abs(lhs), std::abs(rhs));
      INFO("lhs=" << lhs << " rhs=" << rhs << " tol=" << tol);
      REQUIRE(lhs <= rhs + tol);
    };

    const auto require_ge = [&](const Real lhs, const Real rhs) {
      const Real tol = kAbsoluteTol + inequality_tol * max3<Real>(1, std::abs(lhs), std::abs(rhs));
      INFO("lhs=" << lhs << " rhs=" << rhs << " tol=" << tol);
      REQUIRE(lhs + tol >= rhs);
    };

    // Direct cloud formula for active cloud-liquid lanes:
    //   epsc = 2*pi*rho*dv*cdist.
    const auto cloud_expected = [](const RelaxationLane& lane) {
      return 2 * C::Pi * lane.rho * lane.dv * lane.cdist;
    };

    // Analytic rain term with the ventilation term disabled:
    //   epsr = 2*pi*cdistr*rho*dv*f1r*Gamma(mu_r + 2)/lamr.
    const auto analytic_expected = [](const RelaxationLane& lane, const Scalar& f1r) {
      return 2 * C::Pi * lane.cdistr * lane.rho * lane.dv * f1r
             * std::tgamma(lane.mu_r + 2) / lane.lamr;
    };

    // Same analytic rain term, but using the gamma recurrence
    //   Gamma(mu_r + 2) = (mu_r + 1)*Gamma(mu_r + 1).
    //
    // This gives a partly independent check of the implementation formula.
    const auto analytic_recurrence_expected = [](const RelaxationLane& lane, const Scalar& f1r) {
      return 2 * C::Pi * lane.cdistr * lane.rho * lane.dv * f1r
             * (lane.mu_r + 1) * std::tgamma(lane.mu_r + 1) / lane.lamr;
    };

    SECTION("masking_and_thresholds") {
      // Verify activation logic for rain and cloud relaxation.
      //
      // Rain is active only where:
      //   context && qr_incld >= QSMALL.
      //
      // Cloud liquid is active only where:
      //   context && qc_incld >= QSMALL.
      //
      // Important implementation detail:
      // epsr is reset to zero for all lanes, even outside context. epsc is
      // only written inside context, so outside-context epsc should retain
      // its input value.
      auto lanes = make_active_lanes();
      lanes[0].context = false;
      lanes[0].qr_incld = 10 * C::QSMALL;
      lanes[0].qc_incld = 10 * C::QSMALL;
      lanes[0].epsr_in = 777.0;
      lanes[0].epsc_in = -12345.0;

      lanes[1].context = true;
      lanes[1].qr_incld = 0.999 * C::QSMALL;
      lanes[1].qc_incld = 0.999 * C::QSMALL;

      lanes[2].context = true;
      lanes[2].qr_incld = C::QSMALL;
      lanes[2].qc_incld = C::QSMALL;

      lanes[3].context = true;
      lanes[3].qr_incld = 10 * C::QSMALL;
      lanes[3].qc_incld = 0.1 * C::QSMALL;

      run_kernel(revap_table_vals, C::f1r, C::f2r, lanes);

      require_zero(lanes[0].epsr_out);
      require_near(lanes[0].epsc_out, lanes[0].epsc_in);

      require_zero(lanes[1].epsr_out);
      require_zero(lanes[1].epsc_out);

      require_nonnegative(lanes[2].epsr_out);
      require_finite(lanes[2].epsr_out);
      require_near(lanes[2].epsc_out, cloud_expected(lanes[2]));

      require_nonnegative(lanes[3].epsr_out);
      require_finite(lanes[3].epsr_out);
      require_zero(lanes[3].epsc_out);
    }

    SECTION("cloud_identity_and_scaling") {
      // With cloud active and rain disabled, epsc should exactly satisfy:
      //
      //   epsc = 2*pi*rho*dv*cdist.
      //
      // Therefore doubling rho, dv, or cdist should double epsc. The ratio
      // epsc/(rho*dv*cdist) should equal 2*pi.
      auto base = make_active_lanes();
      for (auto& lane : base) {
        lane.qr_incld = 0;
        lane.qc_incld = 5 * C::QSMALL;
      }
      run_kernel(revap_table_vals, C::f1r, C::f2r, base);

      auto rho_scaled = base;
      auto dv_scaled = base;
      auto cdist_scaled = base;
      for (auto& lane : rho_scaled) lane.rho *= 2;
      for (auto& lane : dv_scaled) lane.dv *= 2;
      for (auto& lane : cdist_scaled) lane.cdist *= 2;

      run_kernel(revap_table_vals, C::f1r, C::f2r, rho_scaled);
      run_kernel(revap_table_vals, C::f1r, C::f2r, dv_scaled);
      run_kernel(revap_table_vals, C::f1r, C::f2r, cdist_scaled);

      for (Int s = 0; s < Pack::n; ++s) {
        require_near(base[s].epsc_out, cloud_expected(base[s]));
        require_near(rho_scaled[s].epsc_out, 2 * base[s].epsc_out);
        require_near(dv_scaled[s].epsc_out, 2 * base[s].epsc_out);
        require_near(cdist_scaled[s].epsc_out, 2 * base[s].epsc_out);
        require_near(base[s].epsc_out / (base[s].rho * base[s].dv * base[s].cdist),
                     2 * C::Pi);
      }
    }

    SECTION("rain_additivity") {
      // The rain formula is a sum of two independent terms:
      //
      //   analytic term    proportional to f1r
      //   ventilation term proportional to f2r
      //
      // Therefore:
      //
      //   epsr(f1r, f2r) = epsr(f1r, 0) + epsr(0, f2r).
      //
      // epsc does not depend on f1r or f2r, so it should be unchanged across
      // the three calls.
      auto full = make_active_lanes();
      auto analytic_only = full;
      auto vent_only = full;

      run_kernel(revap_table_vals, C::f1r, C::f2r, full);
      run_kernel(revap_table_vals, C::f1r, 0, analytic_only);
      run_kernel(revap_table_vals, 0, C::f2r, vent_only);

      for (Int s = 0; s < Pack::n; ++s) {
        require_near(full[s].epsr_out,
                     analytic_only[s].epsr_out + vent_only[s].epsr_out);
        require_near(full[s].epsc_out, analytic_only[s].epsc_out);
        require_near(full[s].epsc_out, vent_only[s].epsc_out);
      }
    }

    SECTION("rain_prefactor_scaling") {
      // In the full mixed rain formula, dv and cdistr multiply both the
      // analytic and ventilation terms as shared prefactors.
      //
      // Therefore:
      //
      //   epsr(2*dv)     = 2*epsr(dv)
      //   epsr(2*cdistr) = 2*epsr(cdistr).
      auto base = make_active_lanes();
      auto dv_scaled = base;
      auto cdistr_scaled = base;

      for (auto& lane : dv_scaled) lane.dv *= 2;
      for (auto& lane : cdistr_scaled) lane.cdistr *= 2;

      run_kernel(revap_table_vals, C::f1r, C::f2r, base);
      run_kernel(revap_table_vals, C::f1r, C::f2r, dv_scaled);
      run_kernel(revap_table_vals, C::f1r, C::f2r, cdistr_scaled);

      for (Int s = 0; s < Pack::n; ++s) {
        require_near(dv_scaled[s].epsr_out, 2 * base[s].epsr_out);
        require_near(cdistr_scaled[s].epsr_out, 2 * base[s].epsr_out);
      }
    }

    SECTION("analytic_only") {
      // Disable the ventilation term by setting f2r = 0.
      //
      // The remaining formula is:
      //
      //   epsr = 2*pi*cdistr*rho*dv*f1r*Gamma(mu_r + 2)/lamr.
      //
      // This implies exact scaling with rho and f1r, and inverse scaling
      // with lamr. The test also checks the equivalent gamma-recurrence form.
      auto base = make_active_lanes();
      auto rho_scaled = base;
      auto f1_scaled = base;
      auto lamr_scaled = base;

      for (auto& lane : rho_scaled) lane.rho *= 2;
      for (auto& lane : lamr_scaled) lane.lamr *= 2;

      run_kernel(revap_table_vals, C::f1r, 0, base);
      run_kernel(revap_table_vals, C::f1r, 0, rho_scaled);
      run_kernel(revap_table_vals, 2 * C::f1r, 0, f1_scaled);
      run_kernel(revap_table_vals, C::f1r, 0, lamr_scaled);

      for (Int s = 0; s < Pack::n; ++s) {
        require_near(base[s].epsr_out, analytic_expected(base[s], C::f1r));
        require_near(base[s].epsr_out, analytic_recurrence_expected(base[s], C::f1r));
        require_near(rho_scaled[s].epsr_out, 2 * base[s].epsr_out);
        require_near(f1_scaled[s].epsr_out, 2 * base[s].epsr_out);
        require_near(lamr_scaled[s].epsr_out, 0.5 * base[s].epsr_out);
      }
    }

    SECTION("ventilation_only") {
      // Disable the analytic term by setting f1r = 0.
      //
      // With mu_r and lamr held fixed, the lookup-table value is unchanged.
      // The remaining ventilation term scales as:
      //
      //   epsr proportional to rho^(3/2) * mu^(-1/2) * sc^(1/3) * f2r.
      //
      // Therefore:
      //
      //   rho -> 4*rho gives epsr -> 8*epsr
      //   mu  -> 4*mu  gives epsr -> 0.5*epsr
      //   sc  -> 8*sc  gives epsr -> 2*epsr
      //   f2r -> 2*f2r gives epsr -> 2*epsr.
      auto base = make_active_lanes();
      auto rho_scaled = base;
      auto mu_scaled = base;
      auto sc_scaled = base;
      auto f2_scaled = base;

      for (auto& lane : rho_scaled) lane.rho *= 4;
      for (auto& lane : mu_scaled) lane.mu *= 4;
      for (auto& lane : sc_scaled) lane.sc *= 8;

      run_kernel(revap_table_vals, 0, C::f2r, base);
      run_kernel(revap_table_vals, 0, C::f2r, rho_scaled);
      run_kernel(revap_table_vals, 0, C::f2r, mu_scaled);
      run_kernel(revap_table_vals, 0, C::f2r, sc_scaled);
      run_kernel(revap_table_vals, 0, 2 * C::f2r, f2_scaled);

      for (Int s = 0; s < Pack::n; ++s) {
        require_near(rho_scaled[s].epsr_out, 8 * base[s].epsr_out);
        require_near(mu_scaled[s].epsr_out, 0.5 * base[s].epsr_out);
        require_near(sc_scaled[s].epsr_out, 2 * base[s].epsr_out);
        require_near(f2_scaled[s].epsr_out, 2 * base[s].epsr_out);
      }
    }

    SECTION("mixed_rho_bound") {
      // With both rain terms active, rho appears in two ways:
      //
      //   analytic term    scales as rho
      //   ventilation term scales as rho^(3/2)
      //
      // Therefore doubling rho must increase epsr by at least 2 and by at
      // most 2*sqrt(2), assuming non-negative table values:
      //
      //   2*epsr(rho) <= epsr(2*rho) <= 2*sqrt(2)*epsr(rho).
      auto base = make_active_lanes();
      auto rho_scaled = base;
      for (auto& lane : rho_scaled) lane.rho *= 2;

      run_kernel(revap_table_vals, C::f1r, C::f2r, base);
      run_kernel(revap_table_vals, C::f1r, C::f2r, rho_scaled);

      for (Int s = 0; s < Pack::n; ++s) {
        require_ge(rho_scaled[s].epsr_out, 2 * base[s].epsr_out);
        require_le(rho_scaled[s].epsr_out, 2 * std::sqrt(2.0) * base[s].epsr_out);
      }
    }

    SECTION("zero_coefficients") {
      // Check coefficient locality and shared factors.
      //
      // f1r and f2r only affect rain, so setting both to zero should zero
      // epsr but leave epsc unchanged.
      //
      // cdistr only appears in epsr, so setting cdistr to zero should zero
      // epsr but leave epsc unchanged.
      //
      // cdist only appears in epsc, so setting cdist to zero should zero
      // epsc but leave epsr unchanged.
      //
      // dv appears in both formulas, so setting dv to zero should zero both
      // epsr and epsc.
      auto base = make_active_lanes();
      auto zero_f = make_active_lanes();
      auto zero_cdistr = zero_f;
      auto zero_dv = zero_f;
      auto zero_cdist = zero_f;

      for (auto& lane : zero_cdistr) lane.cdistr = 0;
      for (auto& lane : zero_dv) lane.dv = 0;
      for (auto& lane : zero_cdist) lane.cdist = 0;

      run_kernel(revap_table_vals, C::f1r, C::f2r, base);
      run_kernel(revap_table_vals, 0, 0, zero_f);
      run_kernel(revap_table_vals, C::f1r, C::f2r, zero_cdistr);
      run_kernel(revap_table_vals, C::f1r, C::f2r, zero_dv);
      run_kernel(revap_table_vals, C::f1r, C::f2r, zero_cdist);

      for (Int s = 0; s < Pack::n; ++s) {
        require_zero(zero_f[s].epsr_out);
        require_near(zero_f[s].epsc_out, base[s].epsc_out);

        require_zero(zero_cdistr[s].epsr_out);
        require_near(zero_cdistr[s].epsc_out, base[s].epsc_out);

        require_zero(zero_dv[s].epsr_out);
        require_zero(zero_dv[s].epsc_out);

        require_near(zero_cdist[s].epsr_out, base[s].epsr_out);
        require_zero(zero_cdist[s].epsc_out);
      }
    }

    SECTION("nonnegativity_and_finite") {
      // For physically valid positive inputs and active thresholds, both
      // relaxation timescales should be finite and non-negative. This section
      // samples qr_incld and qc_incld at and above QSMALL to exercise the
      // active threshold boundary without assuming a particular Pack::n.
      auto lanes = make_active_lanes();
      constexpr std::array<Real, 3> qr_vals = {
        C::QSMALL,
        1.001 * C::QSMALL,
        10 * C::QSMALL
      };
      constexpr std::array<Real, 3> qc_vals = {
        C::QSMALL,
        1.001 * C::QSMALL,
        10 * C::QSMALL
      };
      for (Int s = 0; s < Pack::n; ++s) {
        lanes[s].qr_incld = qr_vals[s % qr_vals.size()];
        lanes[s].qc_incld = qc_vals[s % qc_vals.size()];
      }

      run_kernel(revap_table_vals, C::f1r, C::f2r, lanes);

      for (Int s = 0; s < Pack::n; ++s) {
        require_nonnegative(lanes[s].epsr_out);
        require_nonnegative(lanes[s].epsc_out);
        require_finite(lanes[s].epsr_out);
        require_finite(lanes[s].epsc_out);
      }
    }
  }

  void run_bfb()
  {
    // Preserve the original bit-for-bit regression test. The property tests
    // above check mathematical behavior, while this section ensures the
    // production outputs remain identical to the stored baseline.
    auto engine = Base::get_engine();

    const auto& revap_table_vals = this->lookup_tables.revap_table_vals;

    using KTH = KokkosTypes<HostDevice>;

    // Set up input data.
    constexpr Scalar qsmall = C::QSMALL;
    constexpr Scalar qr_small = 0.9 * qsmall;
    constexpr Scalar qr_not_small = 2.0 * qsmall;
    constexpr Scalar qc_small = 0.9 * qsmall;
    constexpr Scalar qc_not_small = 2.0 * qsmall;
    CalcLiqRelaxationData self[max_pack_size];

    // Alternate rain and cloud thresholds across lanes so the BFB case covers
    // both active and inactive branches for epsr and epsc.
    for (Int i = 0; i < max_pack_size; ++i) {
      self[i].randomize(engine);
      self[i].qr_incld = (i % 2) ? qr_small : qr_not_small;
      self[i].qc_incld = ((i/2) % 2) ? qc_small : qc_not_small;
      self[i].f1r = C::f1r;
      self[i].f2r = C::f2r;
    }

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < max_pack_size; ++i) {
        self[i].read(Base::m_ifile);
      }
    }

    // Sync to device
    KTH::view_1d<CalcLiqRelaxationData> self_host("self_host", max_pack_size);
    view_1d<CalcLiqRelaxationData> self_device("self_device", max_pack_size);
    std::copy(&self[0], &self[0] + max_pack_size, self_host.data());
    Kokkos::deep_copy(self_device, self_host);

    // Run the lookup from a kernel and copy results back to host
    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Pack::n;

      // Init pack inputs
      Pack rho, dv, mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld;

      for (Int s = 0, vs = offset; s < Pack::n; ++s, ++vs) {
        rho[s]      = self_device(vs).rho;
        dv[s]       = self_device(vs).dv;
        mu[s]       = self_device(vs).mu;
        sc[s]       = self_device(vs).sc;
        mu_r[s]     = self_device(vs).mu_r;
        lamr[s]     = self_device(vs).lamr;
        cdistr[s]   = self_device(vs).cdistr;
        cdist[s]    = self_device(vs).cdist;
        qr_incld[s] = self_device(vs).qr_incld;
        qc_incld[s] = self_device(vs).qc_incld;
      }

      Pack epsr{0.0}, epsc{0.0};
      Functions::calc_liq_relaxation_timescale(revap_table_vals, rho, self_device(0).f1r, self_device(0).f2r, dv,
        mu, sc, mu_r, lamr, cdistr, cdist, qr_incld, qc_incld, epsr, epsc);

      for (Int s = 0, vs = offset; s < Pack::n; ++s, ++vs) {
        self_device(vs).epsr = epsr[s];
        self_device(vs).epsc = epsc[s];
      }
    });

    Kokkos::deep_copy(self_host, self_device);

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        REQUIRE(self[s].epsr == self_host(s).epsr);
        REQUIRE(self[s].epsc == self_host(s).epsc);
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int s = 0; s < max_pack_size; ++s) {
        self_host(s).write(Base::m_ofile);
      }
    }
  }

};

}
}
}

namespace {

TEST_CASE("p3_calc_liq_relaxation_timescale", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcLiqRelaxationTimescale;

  T t;
  SECTION("properties") {
    t.run_phys();
  }

  SECTION("bfb") {
    t.run_bfb();
  }
}

}
