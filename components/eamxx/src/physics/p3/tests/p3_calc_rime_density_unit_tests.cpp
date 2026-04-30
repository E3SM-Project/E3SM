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
 * Property tests for calc_rime_density.
 *
 * For active cloud-liquid riming lanes,
 *
 *   vtrmi1 = table_val_qi_fallspd * rhofaci
 *
 * and, when qc_incld >= QSMALL,
 *
 *   Vt_qc    = acn * (mu_c + 4) * (mu_c + 5) / lambda_c^2
 *   D_c      = (mu_c + 4) / lamc
 *   V_impact = |vtrmi1 - Vt_qc|
 *   inv_Tc   = 1 / min(-0.001, T_atm - T_zerodegc)
 *   Ri       = clamp(-0.5e6 * D_c * V_impact * inv_Tc, 1, 12)
 *
 * with rho_qm_cloud given by a piecewise function of Ri.
 *
 * These tests validate the activation mask, exact local identities, targeted
 * Ri regimes, and the monotone/bounded microphysical behavior implied by the
 * closure, while keeping the production kernel as the source of truth.
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
struct UnitWrap::UnitTest<D>::TestCalcRimeDensity : public UnitWrap::UnitTest<D>::Base {

  struct RimeLane {
    Real T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c;
    Real qc_incld, qc2qi_collect_tend;
    bool context;
    Real vtrmi1_in, rho_qm_cloud_in;
    Real vtrmi1_out, rho_qm_cloud_out;
  };

  std::array<RimeLane, max_pack_size> make_active_lanes() const
  {
    std::array<RimeLane, max_pack_size> lanes;
    for (Int s = 0; s < max_pack_size; ++s) {
      lanes[s] = {
        C::T_zerodegc.value - (5.0 + 0.5 * s),
        0.65 + 0.03 * s,
        0.8 + 0.08 * s,
        0.35 + 0.04 * s,
        2.8 + 0.15 * s,
        1.1 + 0.08 * s,
        5.0 * C::QSMALL,
        4.0 * C::QSMALL,
        true,
        0.0,
        400.0,
        0.0,
        0.0
      };
    }
    return lanes;
  }

  static Scalar cloud_droplet_fall_speed(const Scalar acn, const Scalar lamc,
                                         const Scalar mu_c)
  {
    return acn * (mu_c + 4) * (mu_c + 5) / (lamc * lamc);
  }

  static Scalar density_from_ri(const Scalar ri_raw)
  {
    const Scalar ri = std::max(1.0, std::min(12.0, ri_raw));
    return ri <= 8.0
      ? 1000.0 * (0.051 + 0.114 * ri - 0.0055 * ri * ri)
      : 611.0 + 72.25 * (ri - 8.0);
  }

  static Scalar ri_from_state(const RimeLane& lane)
  {
    const Scalar vtrmi1 = lane.table_val_qi_fallspd * lane.rhofaci;
    const Scalar vt_qc = cloud_droplet_fall_speed(lane.acn, lane.lamc, lane.mu_c);
    const Scalar d_c = (lane.mu_c + 4) / lane.lamc;
    const Scalar v_impact = std::abs(vtrmi1 - vt_qc);
    const Scalar denom = std::min(-0.001, lane.T_atm - C::T_zerodegc.value);
    const Scalar ri_raw = -0.5e6 * d_c * v_impact / denom;
    return std::max(1.0, std::min(12.0, ri_raw));
  }

  static Scalar active_density_expected(const RimeLane& lane)
  {
    return density_from_ri(ri_from_state(lane));
  }

  static Scalar table_val_for_target_ri(const Scalar ri_raw, const Scalar T_atm,
                                        const Scalar rhofaci, const Scalar acn,
                                        const Scalar lamc, const Scalar mu_c)
  {
    const Scalar delta_t = C::T_zerodegc.value - T_atm;
    const Scalar effective_delta_t = std::max<Scalar>(0.001, delta_t);
    const Scalar d_c = (mu_c + 4) / lamc;
    const Scalar v_impact = ri_raw * effective_delta_t / (5.0e5 * d_c);
    const Scalar vt_qc = cloud_droplet_fall_speed(acn, lamc, mu_c);
    const Scalar vtrmi1 = vt_qc + v_impact;
    return vtrmi1 / rhofaci;
  }

  void run_kernel(std::array<RimeLane, max_pack_size>& lanes) const
  {
    using KTH = KokkosTypes<HostDevice>;

    KTH::view_1d<RimeLane> lanes_host("lanes_host", max_pack_size);
    view_1d<RimeLane> lanes_device("lanes_device", max_pack_size);
    std::copy(lanes.begin(), lanes.end(), lanes_host.data());
    Kokkos::deep_copy(lanes_device, lanes_host);

    Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
      const Int offset = i * Pack::n;
      Pack T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c;
      Pack qc_incld, qc2qi_collect_tend, vtrmi1, rho_qm_cloud;
      Mask context;

      for (Int s = 0; s < Pack::n; ++s) {
        const Int vs = offset + s;
        T_atm[s] = lanes_device(vs).T_atm;
        rhofaci[s] = lanes_device(vs).rhofaci;
        table_val_qi_fallspd[s] = lanes_device(vs).table_val_qi_fallspd;
        acn[s] = lanes_device(vs).acn;
        lamc[s] = lanes_device(vs).lamc;
        mu_c[s] = lanes_device(vs).mu_c;
        qc_incld[s] = lanes_device(vs).qc_incld;
        qc2qi_collect_tend[s] = lanes_device(vs).qc2qi_collect_tend;
        vtrmi1[s] = lanes_device(vs).vtrmi1_in;
        rho_qm_cloud[s] = lanes_device(vs).rho_qm_cloud_in;
        context.set(s, lanes_device(vs).context);
      }

      Functions::calc_rime_density(T_atm, rhofaci, table_val_qi_fallspd, acn,
                                   lamc, mu_c, qc_incld, qc2qi_collect_tend,
                                   vtrmi1, rho_qm_cloud, context);

      for (Int s = 0; s < Pack::n; ++s) {
        const Int vs = offset + s;
        lanes_device(vs).vtrmi1_out = vtrmi1[s];
        lanes_device(vs).rho_qm_cloud_out = rho_qm_cloud[s];
      }
    });

    Kokkos::deep_copy(lanes_host, lanes_device);
    std::copy(lanes_host.data(), lanes_host.data() + max_pack_size, lanes.begin());
  }

void run_phys()
{
  const Real identity_tol = 100 * std::numeric_limits<Real>::epsilon();
  const Real inequality_tol = 100 * std::numeric_limits<Real>::epsilon();

  const auto near_tol = [&](const Real actual, const Real expected) {
    return kAbsoluteTol + identity_tol * max3<Real>(1, std::abs(actual), std::abs(expected));
  };

  const auto require_near = [&](const Real actual, const Real expected) {
    INFO("actual=" << actual << " expected=" << expected
         << " tol=" << near_tol(actual, expected));
    REQUIRE(std::abs(actual - expected) <= near_tol(actual, expected));
  };

  const auto require_density_fit = [&](const Real actual, const Real expected) {
    const Real density_fit_tol = std::max<Real>(
      1e-6,
      kAbsoluteTol + 200 * std::numeric_limits<Real>::epsilon()
        * max3<Real>(1, std::abs(actual), std::abs(expected)));
    INFO("actual=" << actual << " expected=" << expected
         << " abs_tol=" << density_fit_tol);
    REQUIRE(std::abs(actual - expected) <= density_fit_tol);
  };

  const auto require_zero = [&](const Real actual) {
    INFO("actual=" << actual << " abs_tol=" << kAbsoluteTol);
    REQUIRE(std::abs(actual) <= kAbsoluteTol);
  };

  const auto require_finite = [&](const Real actual) {
    INFO("actual=" << actual);
    REQUIRE(std::isfinite(actual));
  };

  const auto require_ge = [&](const Real lhs, const Real rhs) {
    const Real tol = kAbsoluteTol + inequality_tol * max3<Real>(1, std::abs(lhs), std::abs(rhs));
    INFO("lhs=" << lhs << " rhs=" << rhs << " tol=" << tol);
    REQUIRE(lhs + tol >= rhs);
  };

  const auto require_le = [&](const Real lhs, const Real rhs) {
    const Real tol = kAbsoluteTol + inequality_tol * max3<Real>(1, std::abs(lhs), std::abs(rhs));
    INFO("lhs=" << lhs << " rhs=" << rhs << " tol=" << tol);
    REQUIRE(lhs <= rhs + tol);
  };

  SECTION("activation_and_defaults") {
    // The active mask requires qc2qi_collect_tend >= QSMALL, T_atm below
    // freezing, and context = true. qc_incld only controls whether the
    // density formula is used or the 400 kg m^-3 default is retained.
    auto lanes = make_active_lanes();

    lanes[0].qc2qi_collect_tend = 0.9 * C::QSMALL;
    lanes[1].qc2qi_collect_tend = C::QSMALL;
    lanes[2].qc_incld = 0.9 * C::QSMALL;
    lanes[3].qc_incld = C::QSMALL;
    lanes[4].T_atm = C::T_zerodegc.value;
    lanes[5].T_atm = C::T_zerodegc.value - 1.0;

    run_kernel(lanes);

    require_zero(lanes[0].vtrmi1_out);
    require_near(lanes[0].rho_qm_cloud_out, 400.0);

    require_near(lanes[1].vtrmi1_out, lanes[1].table_val_qi_fallspd * lanes[1].rhofaci);
    require_near(lanes[1].rho_qm_cloud_out, active_density_expected(lanes[1]));

    require_near(lanes[2].vtrmi1_out, lanes[2].table_val_qi_fallspd * lanes[2].rhofaci);
    require_near(lanes[2].rho_qm_cloud_out, 400.0);

    require_near(lanes[3].vtrmi1_out, lanes[3].table_val_qi_fallspd * lanes[3].rhofaci);
    require_near(lanes[3].rho_qm_cloud_out, active_density_expected(lanes[3]));

    require_zero(lanes[4].vtrmi1_out);
    require_near(lanes[4].rho_qm_cloud_out, 400.0);

    require_near(lanes[5].vtrmi1_out, lanes[5].table_val_qi_fallspd * lanes[5].rhofaci);
    require_near(lanes[5].rho_qm_cloud_out, active_density_expected(lanes[5]));
  }

  SECTION("velocity_identity") {
    // In every active collection lane, the mean ice fall speed used in riming
    // is exactly the lookup-table fall speed scaled by rhofaci.
    auto lanes = make_active_lanes();
    run_kernel(lanes);

    for (Int s = 0; s < max_pack_size; ++s) {
      require_near(lanes[s].vtrmi1_out, lanes[s].table_val_qi_fallspd * lanes[s].rhofaci);
    }
  }

  SECTION("targeted_ri_density_regimes") {
    // Construct lanes by solving backward for table_val_qi_fallspd from a
    // target unclamped Ri value so we can hit the low clamp, branch point, and
    // wet-growth upper clamp deterministically.
    auto lanes = make_active_lanes();
    constexpr std::array<Real, 7> ri_raw = {0.5, 1.0, 4.0, 8.0, 10.0, 12.0, 20.0};
    constexpr std::array<Real, 7> rho_expected = {159.5, 159.5, 419.0, 611.0, 755.5, 900.0, 900.0};
    const Real T_atm = C::T_zerodegc.value - 1.0;
    const Real rhofaci = 0.8;
    const Real acn = 0.5;
    const Real lamc = 3.5;
    const Real mu_c = 1.2;

    for (Int i = 0; i < static_cast<Int>(ri_raw.size()); ++i) {
      lanes[i].T_atm = T_atm;
      lanes[i].rhofaci = rhofaci;
      lanes[i].acn = acn;
      lanes[i].lamc = lamc;
      lanes[i].mu_c = mu_c;
      lanes[i].table_val_qi_fallspd = table_val_for_target_ri(ri_raw[i], T_atm, rhofaci, acn, lamc, mu_c);
    }

    run_kernel(lanes);

    for (Int i = 0; i < static_cast<Int>(ri_raw.size()); ++i) {
      require_density_fit(lanes[i].rho_qm_cloud_out, rho_expected[i]);
      require_near(lanes[i].vtrmi1_out, lanes[i].table_val_qi_fallspd * lanes[i].rhofaci);
    }
  }

  SECTION("rime_density_increases_with_impact_speed") {
    // At fixed droplet properties and temperature, Ri is proportional to
    // impact speed, so rho_qm_cloud should be nondecreasing with Ri and plateau
    // at the low and high clamp limits.
    auto lanes = make_active_lanes();
    constexpr std::array<Real, 8> ri_raw = {0.5, 1.0, 2.0, 4.0, 8.0, 10.0, 12.0, 20.0};
    const Real T_atm = C::T_zerodegc.value - 1.0;
    const Real rhofaci = 0.8;
    const Real acn = 0.5;
    const Real lamc = 3.5;
    const Real mu_c = 1.2;

    for (Int i = 0; i < static_cast<Int>(ri_raw.size()); ++i) {
      lanes[i].T_atm = T_atm;
      lanes[i].rhofaci = rhofaci;
      lanes[i].acn = acn;
      lanes[i].lamc = lamc;
      lanes[i].mu_c = mu_c;
      lanes[i].table_val_qi_fallspd = table_val_for_target_ri(ri_raw[i], T_atm, rhofaci, acn, lamc, mu_c);
    }

    run_kernel(lanes);

    for (Int i = 1; i < static_cast<Int>(ri_raw.size()); ++i) {
      require_ge(lanes[i].rho_qm_cloud_out, lanes[i-1].rho_qm_cloud_out);
    }
    require_density_fit(lanes[0].rho_qm_cloud_out, lanes[1].rho_qm_cloud_out);
    require_density_fit(lanes[6].rho_qm_cloud_out, lanes[7].rho_qm_cloud_out);
  }

  SECTION("rime_density_increases_toward_freezing") {
    // For fixed droplet size and impact speed, the closure increases Ri as the
    // temperature approaches freezing from below, so rho_qm_cloud should be
    // nondecreasing along that sequence.
    auto lanes = make_active_lanes();
    constexpr std::array<Real, 4> temps = {
      C::T_zerodegc.value - 20.0,
      C::T_zerodegc.value - 10.0,
      C::T_zerodegc.value - 1.0,
      C::T_zerodegc.value - 0.001
    };
    const Real rhofaci = 0.8;
    const Real acn = 0.5;
    const Real lamc = 3.5;
    const Real mu_c = 1.2;
    const Real table_val = table_val_for_target_ri(4.0, C::T_zerodegc.value - 1.0,
                                                   rhofaci, acn, lamc, mu_c);

    for (Int i = 0; i < static_cast<Int>(temps.size()); ++i) {
      lanes[i].T_atm = temps[i];
      lanes[i].rhofaci = rhofaci;
      lanes[i].acn = acn;
      lanes[i].lamc = lamc;
      lanes[i].mu_c = mu_c;
      lanes[i].table_val_qi_fallspd = table_val;
    }

    run_kernel(lanes);

    for (Int i = 1; i < static_cast<Int>(temps.size()); ++i) {
      require_ge(lanes[i].rho_qm_cloud_out, lanes[i-1].rho_qm_cloud_out);
    }
  }

  SECTION("impact_speed_symmetry") {
    // The closure depends on V_impact through an absolute value, so equal
    // positive perturbations above and below Vt_qc should produce the same
    // rime density.
    auto lanes = make_active_lanes();
    const Real T_atm = C::T_zerodegc.value - 5.0;
    const Real rhofaci = 0.8;
    const Real acn = 0.5;
    const Real lamc = 3.5;
    const Real mu_c = 1.2;
    const Real vt_qc = cloud_droplet_fall_speed(acn, lamc, mu_c);
    const Real delta = 0.1 * vt_qc;
    const Real vtrmi1_hi = vt_qc + delta;
    const Real vtrmi1_lo = vt_qc - delta;

    lanes[0].T_atm = T_atm;
    lanes[0].rhofaci = rhofaci;
    lanes[0].acn = acn;
    lanes[0].lamc = lamc;
    lanes[0].mu_c = mu_c;
    lanes[0].table_val_qi_fallspd = vtrmi1_hi / rhofaci;

    lanes[1].T_atm = T_atm;
    lanes[1].rhofaci = rhofaci;
    lanes[1].acn = acn;
    lanes[1].lamc = lamc;
    lanes[1].mu_c = mu_c;
    lanes[1].table_val_qi_fallspd = vtrmi1_lo / rhofaci;

    run_kernel(lanes);

    require_near(lanes[0].rho_qm_cloud_out, lanes[1].rho_qm_cloud_out);
  }

  SECTION("rime_density_bounds") {
    // Active cloud-riming lanes with qc_incld >= QSMALL should stay within the
    // closure bounds, while inactive context-true lanes should return the 400
    // kg m^-3 default.
    auto lanes = make_active_lanes();
    constexpr std::array<Real, 7> ri_raw = {0.5, 1.0, 4.0, 8.0, 10.0, 12.0, 20.0};
    const Real T_atm = C::T_zerodegc.value - 1.0;
    const Real rhofaci = 0.8;
    const Real acn = 0.5;
    const Real lamc = 3.5;
    const Real mu_c = 1.2;

    for (Int i = 0; i < static_cast<Int>(ri_raw.size()); ++i) {
      lanes[i].T_atm = T_atm;
      lanes[i].rhofaci = rhofaci;
      lanes[i].acn = acn;
      lanes[i].lamc = lamc;
      lanes[i].mu_c = mu_c;
      lanes[i].table_val_qi_fallspd = table_val_for_target_ri(ri_raw[i], T_atm, rhofaci, acn, lamc, mu_c);
    }
    lanes[7].qc2qi_collect_tend = 0.9 * C::QSMALL;

    run_kernel(lanes);

    for (Int i = 0; i < static_cast<Int>(ri_raw.size()); ++i) {
      require_ge(lanes[i].rho_qm_cloud_out, 159.5);
      require_le(lanes[i].rho_qm_cloud_out, 900.0);
    }
    require_zero(lanes[7].vtrmi1_out);
    require_near(lanes[7].rho_qm_cloud_out, 400.0);
  }

  SECTION("temperature_floor_prevents_singularity") {
    // The denominator in inv_Tc is floored at 0.001 K below freezing, so the
    // closure should remain finite and identical for temperatures at or above
    // that floor while remaining inactive at and above freezing.
    auto lanes = make_active_lanes();
    const Real rhofaci = 0.8;
    const Real acn = 0.5;
    const Real lamc = 3.5;
    const Real mu_c = 1.2;
    const Real table_val = table_val_for_target_ri(10.0, C::T_zerodegc.value - 0.001,
                                                   rhofaci, acn, lamc, mu_c);

    lanes[0].T_atm = C::T_zerodegc.value - 0.001;
    lanes[1].T_atm = C::T_zerodegc.value - 0.0001;
    lanes[2].T_atm = C::T_zerodegc.value;
    lanes[3].T_atm = C::T_zerodegc.value + 1.0;
    for (Int i = 0; i < 4; ++i) {
      lanes[i].rhofaci = rhofaci;
      lanes[i].acn = acn;
      lanes[i].lamc = lamc;
      lanes[i].mu_c = mu_c;
      lanes[i].table_val_qi_fallspd = table_val;
    }

    run_kernel(lanes);

    require_finite(lanes[0].rho_qm_cloud_out);
    require_finite(lanes[1].rho_qm_cloud_out);
    require_ge(lanes[0].rho_qm_cloud_out, 159.5);
    require_le(lanes[0].rho_qm_cloud_out, 900.0);
    require_ge(lanes[1].rho_qm_cloud_out, 159.5);
    require_le(lanes[1].rho_qm_cloud_out, 900.0);
    require_near(lanes[0].rho_qm_cloud_out, lanes[1].rho_qm_cloud_out);
    require_near(lanes[0].vtrmi1_out, lanes[1].vtrmi1_out);

    require_zero(lanes[2].vtrmi1_out);
    require_near(lanes[2].rho_qm_cloud_out, 400.0);
    require_zero(lanes[3].vtrmi1_out);
    require_near(lanes[3].rho_qm_cloud_out, 400.0);
  }

  SECTION("context_mask_preserves_inactive_lanes") {
    // When context is false, the kernel should not overwrite either output, so
    // sentinel inputs let us verify that masked lanes are preserved.
    auto lanes = make_active_lanes();
    for (auto& lane : lanes) {
      lane.context = false;
      lane.vtrmi1_in = -999.0;
      lane.rho_qm_cloud_in = -777.0;
    }

    run_kernel(lanes);

    for (Int s = 0; s < max_pack_size; ++s) {
      require_near(lanes[s].vtrmi1_out, lanes[s].vtrmi1_in);
      require_near(lanes[s].rho_qm_cloud_out, lanes[s].rho_qm_cloud_in);
    }
  }
}

void run_bfb()
{
  // This is the threshold for whether the qc cloud mixing ratio is large
  // enough to affect the rime density.
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar qc_incld_small = 0.9 * qsmall;
  constexpr Scalar qc_incld_not_small = 2.0 * qsmall;

  // Same goes for qc2qi_collect_tend.
  constexpr Scalar qc2qi_collect_tend_small = 0.9 * qsmall;
  constexpr Scalar qc2qi_collect_tend_not_small = 2.0 * qsmall;

  // We need to test the calculation under freezing and not-freezing
  // conditions.
  constexpr Scalar t_freezing = 0.9 * C::T_rainfrz.value,
                   t_not_freezing = 2.0 * C::T_rainfrz.value;

  // The deterministic property tests above target Ri directly. This BFB section
  // preserves the historical regression cases and baseline layout.
  constexpr Scalar rhofaci1 = 0.25, rhofaci2 = 0.5;
  constexpr Scalar table_val_qi_fallspd1 = 0.125, table_val_qi_fallspd2 = 0.875;
  constexpr Scalar acn1 = 0.1, acn2 = 0.6;
  constexpr Scalar lamc1 = 0.4, lamc2 = 0.8;
  constexpr Scalar mu_c1 = 0.2, mu_c2 = 0.7;

  CalcRimeDensityData calc_rime_density_data[max_pack_size] = {
    // T_atm, rhofaci, table_val_qi_fallspd1, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend
    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_small},

    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_not_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_not_small},

    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_small},

    {t_not_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_not_small},
    {t_not_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_not_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_small, qc2qi_collect_tend_not_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_small, qc2qi_collect_tend_not_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_small},

    {t_freezing, rhofaci1, table_val_qi_fallspd1, acn1, lamc1, mu_c1, qc_incld_not_small, qc2qi_collect_tend_not_small},
    {t_freezing, rhofaci2, table_val_qi_fallspd2, acn2, lamc2, mu_c2, qc_incld_not_small, qc2qi_collect_tend_not_small},
  };

  // Sync to device
  view_1d<CalcRimeDensityData> device_data("calc_rime_density", max_pack_size);
  const auto host_data = Kokkos::create_mirror_view(device_data);
  std::copy(&calc_rime_density_data[0], &calc_rime_density_data[0] + max_pack_size,
            host_data.data());
  Kokkos::deep_copy(device_data, host_data);

  // Read baseline data
  if (this->m_baseline_action == COMPARE) {
    for (Int i = 0; i < max_pack_size; ++i) {
      calc_rime_density_data[i].read(Base::m_ifile);
    }
  }

  // Run the lookup from a kernel and copy results back to host
  Kokkos::parallel_for(num_test_itrs, KOKKOS_LAMBDA(const Int& i) {
    const Int offset = i * Pack::n;

    // Init pack inputs
    Pack T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c, qc_incld, qc2qi_collect_tend;
    for (Int s = 0, vs = offset; s < Pack::n; ++s, ++vs) {
      T_atm[s]                = device_data(vs).T_atm;
      rhofaci[s]              = device_data(vs).rhofaci;
      table_val_qi_fallspd[s] = device_data(vs).table_val_qi_fallspd;
      acn[s]                  = device_data(vs).acn;
      lamc[s]                 = device_data(vs).lamc;
      mu_c[s]                 = device_data(vs).mu_c;
      qc_incld[s]             = device_data(vs).qc_incld;
      qc2qi_collect_tend[s]   = device_data(vs).qc2qi_collect_tend;
    }

    Pack vtrmi1{0.0};
    Pack rho_qm_cloud{0.0};

    Functions::calc_rime_density(T_atm, rhofaci, table_val_qi_fallspd, acn, lamc, mu_c,
                                 qc_incld, qc2qi_collect_tend, vtrmi1, rho_qm_cloud);

    // Copy results back into views
    for (Int s = 0, vs = offset; s < Pack::n; ++s, ++vs) {
      device_data(vs).vtrmi1  = vtrmi1[s];
      device_data(vs).rho_qm_cloud  = rho_qm_cloud[s];
    }
  });

  // Sync back to host.
  Kokkos::deep_copy(host_data, device_data);

  // Validate results.
  if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      REQUIRE(calc_rime_density_data[s].vtrmi1 == host_data[s].vtrmi1);
      REQUIRE(calc_rime_density_data[s].rho_qm_cloud == host_data[s].rho_qm_cloud);
    }
  }
  else if (this->m_baseline_action == GENERATE) {
    for (Int s = 0; s < max_pack_size; ++s) {
      host_data(s).write(Base::m_ofile);
    }
  }
}

};

}
}
}

namespace {

TEST_CASE("p3_calc_rime_density", "[p3_functions]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcRimeDensity;

  T t;
  SECTION("properties") {
    t.run_phys();
  }

  SECTION("bfb") {
    t.run_bfb();
  }
}

} // namespace
