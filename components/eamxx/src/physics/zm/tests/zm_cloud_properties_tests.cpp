#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmCloudProperties : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmCloudPropertiesData baseline_data[] = {
      //                    pcols, ncol, pver, pverp, msg, limcnv
      ZmCloudPropertiesData(    4,    4,   72,    73,  10,  40),
      ZmCloudPropertiesData(    4,    4,   72,    73,  20,  50),
      ZmCloudPropertiesData(    4,    4,  128,   129,  30,  60),
      ZmCloudPropertiesData(    4,    4,  128,   129,  40,  70),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmCloudPropertiesData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmCloudPropertiesData test_data[] = {
      ZmCloudPropertiesData(baseline_data[0]),
      ZmCloudPropertiesData(baseline_data[1]),
      ZmCloudPropertiesData(baseline_data[2]),
      ZmCloudPropertiesData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      zm_cloud_properties(d);
    }

    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmCloudPropertiesData& d_baseline = baseline_data[i];
        ZmCloudPropertiesData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.mflx_up));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.entr_up));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.detr_up));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.mflx_dn));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.entr_dn));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.mflx_net));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.s_upd));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.q_upd));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.ql));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.s_dnd));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.q_dnd));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.qst));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.cu));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.evp));
        REQUIRE(d_baseline.total(d_baseline.mflx_up) == d_test.total(d_test.rprd));
        for (Int k = 0; k < d_baseline.total(d_baseline.mflx_up); ++k) {
          REQUIRE(d_baseline.mflx_up[k] == Approx(d_test.mflx_up[k]).margin(margin));
          REQUIRE(d_baseline.entr_up[k] == Approx(d_test.entr_up[k]).margin(margin));
          REQUIRE(d_baseline.detr_up[k] == Approx(d_test.detr_up[k]).margin(margin));
          REQUIRE(d_baseline.mflx_dn[k] == Approx(d_test.mflx_dn[k]).margin(margin));
          REQUIRE(d_baseline.entr_dn[k] == Approx(d_test.entr_dn[k]).margin(margin));
          REQUIRE(d_baseline.mflx_net[k] == Approx(d_test.mflx_net[k]).margin(margin));
          REQUIRE(d_baseline.s_upd[k] == Approx(d_test.s_upd[k]).margin(margin));
          REQUIRE(d_baseline.q_upd[k] == Approx(d_test.q_upd[k]).margin(margin));
          REQUIRE(d_baseline.ql[k] == Approx(d_test.ql[k]).margin(margin));
          REQUIRE(d_baseline.s_dnd[k] == Approx(d_test.s_dnd[k]).margin(margin));
          REQUIRE(d_baseline.q_dnd[k] == Approx(d_test.q_dnd[k]).margin(margin));
          REQUIRE(d_baseline.qst[k] == Approx(d_test.qst[k]).margin(margin));
          REQUIRE(d_baseline.cu[k] == Approx(d_test.cu[k]).margin(margin));
          REQUIRE(d_baseline.evp[k] == Approx(d_test.evp[k]).margin(margin));
          REQUIRE(d_baseline.rprd[k] == Approx(d_test.rprd[k]).margin(margin));
        }
        REQUIRE(d_baseline.total(d_baseline.pflx) == d_test.total(d_test.pflx));
        for (Int k = 0; k < d_baseline.total(d_baseline.pflx); ++k) {
          REQUIRE(d_baseline.pflx[k] == Approx(d_test.pflx[k]).margin(margin));
        }
        REQUIRE(d_baseline.total(d_baseline.jt) == d_test.total(d_test.jt));
        REQUIRE(d_baseline.total(d_baseline.jt) == d_test.total(d_test.jlcl));
        REQUIRE(d_baseline.total(d_baseline.jt) == d_test.total(d_test.j0));
        REQUIRE(d_baseline.total(d_baseline.jt) == d_test.total(d_test.jd));
        for (Int k = 0; k < d_baseline.total(d_baseline.jt); ++k) {
          REQUIRE(d_baseline.jt[k] == d_test.jt[k]);
          REQUIRE(d_baseline.jlcl[k] == d_test.jlcl[k]);
          REQUIRE(d_baseline.j0[k] == d_test.j0[k]);
          REQUIRE(d_baseline.jd[k] == d_test.jd[k]);
        }
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        test_data[i].write(Base::m_ofile);
      }
    }
  } // run_bfb

  // Property test for the turbulent entrainment/detrainment offset
  // (frac_turb_mix). This is deliberately NOT a "bfb" test: it compares a
  // baseline run (frac_turb_mix=0) against a run with a positive offset on
  // identical inputs and checks the invariants the offset must satisfy,
  // rather than comparing against a stored Fortran baseline.
  void run_property()
  {
    auto engine = Base::get_engine();

    // Positive turbulent mixing rate [1/m], comparable in magnitude to the
    // dynamic entrainment parameter dmpdz and large enough to exercise the
    // positivity cap that bounds the added detrainment.
    const Real frac_turb_mix = 1e-3;

    ZmCloudPropertiesData base_data[] = {
      //                    pcols, ncol, pver, pverp, msg, limcnv
      ZmCloudPropertiesData(    4,    4,   72,    73,  10,  40),
      ZmCloudPropertiesData(    4,    4,   72,    73,  20,  50),
      ZmCloudPropertiesData(    4,    4,  128,   129,  30,  60),
      ZmCloudPropertiesData(    4,    4,  128,   129,  40,  70),
    };
    static constexpr Int num_runs = sizeof(base_data) / sizeof(ZmCloudPropertiesData);

    for (auto& d : base_data) {
      d.randomize(engine);
    }

    // Copies run with the turbulent mixing enabled. The copy constructor does
    // not carry frac_turb_mix (it is not a PTD_STD_DEF scalar), so it must be
    // set explicitly here after the copy.
    ZmCloudPropertiesData mix_data[] = {
      ZmCloudPropertiesData(base_data[0]),
      ZmCloudPropertiesData(base_data[1]),
      ZmCloudPropertiesData(base_data[2]),
      ZmCloudPropertiesData(base_data[3]),
    };
    for (auto& d : mix_data) {
      d.frac_turb_mix = frac_turb_mix;
    }

    for (auto& d : base_data) { zm_cloud_properties(d); }
    for (auto& d : mix_data)  { zm_cloud_properties(d); }

    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 10000 : 100);

    // Aggregated across every config/column: used to confirm the offset
    // actually changes the solution somewhere convection is active.
    Real max_entr_diff = 0;
    Real total_active  = 0;

    for (Int i = 0; i < num_runs; ++i) {
      ZmCloudPropertiesData& b = base_data[i];
      ZmCloudPropertiesData& m = mix_data[i];
      const Int pver  = m.pver;
      const Int pverp = m.pverp;
      const Int pcols = m.pcols;

      for (Int c = 0; c < pcols; ++c) {
        for (Int k = 0; k < pver; ++k) {
          const Int idx = c*pver + k;

          // (1) Nonnegative updraft condensate: the positivity cap must keep
          //     the enhanced detrainment from driving liquid water negative.
          REQUIRE(m.ql[idx] >= -margin);

          // Accumulators for the mixing-has-an-effect check below.
          total_active  += std::abs(m.mflx_up[idx]);
          max_entr_diff  = std::max(max_entr_diff,
                                    std::abs(m.entr_up[idx] - b.entr_up[idx]));

          // (2) Equal increments / net mass flux preserved: adding the same
          //     offset to entrainment and detrainment must leave the mass-flux
          //     divergence identity (E-D)*dz == dMflx intact. An unequal
          //     offset would break this.
          if (k < pver - 1) {
            const Int zidx = c*pverp + k;
            const Real dz  = m.z_int[zidx] - m.z_int[zidx+1];
            const Real lhs = (m.entr_up[idx] - m.detr_up[idx]) * dz;
            const Real rhs =  m.mflx_up[idx] - m.mflx_up[idx+1];
            REQUIRE(lhs == Approx(rhs).margin(margin));
          }
        }
      }
    }

    // (3) The offset must actually change the solution somewhere convection is
    //     active - this catches a regression that drops it entirely. Pooled
    //     across all configs/columns so that a single thin/inactive cloud
    //     cannot make the check vacuous.
    if (total_active > 0) {
      REQUIRE(max_entr_diff > 0);
    }
  } // run_property

};

} // namespace unit_test
} // namespace zm
} // namespace scream

namespace {

TEST_CASE("zm_cloud_properties_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmCloudProperties;

  TestStruct t;
  t.run_bfb();
}

TEST_CASE("zm_cloud_properties_turb_mix", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmCloudProperties;

  TestStruct t;
  t.run_property();
}

} // empty namespace
