#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmConvMain : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs. NOTE: We only due single-column tests in order to avoid
    // issues with gather/cross-column dependencies that we aren't going to deal
    // with on the eamxx side.
    ZmConvMainData baseline_data[] = {
      //             pcols, ncol, pver, pverp, time_step, is_first_step, lengath (output, set to zero)
      ZmConvMainData(    1,    1,   72,    73,       2.0,         true, 0),
      ZmConvMainData(    1,    1,   72,    73,       3.0,         false , 0),
      ZmConvMainData(    1,    1,  128,   129,       4.0,         true, 0),
      ZmConvMainData(    1,    1,  128,   129,       5.0,         false , 0),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmConvMainData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmConvMainData test_data[] = {
      ZmConvMainData(baseline_data[0]),
      ZmConvMainData(baseline_data[1]),
      ZmConvMainData(baseline_data[2]),
      ZmConvMainData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      if (this->m_baseline_action == GENERATE) {
        zm_conv_main_f(d);
      }
      else {
        zm_conv_main(d);
      }
    }

    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmConvMainData& d_baseline = baseline_data[i];
        ZmConvMainData& d_test = test_data[i];
        REQUIRE(d_baseline.lengath == d_test.lengath);
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.prec));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.cape));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.dcape));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.dsubcld));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.rliq));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.gather_index));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.msemax_klev));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.jctop));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.jcbot));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.jt));
        for (Int k = 0; k < d_baseline.total(d_baseline.prec); ++k) {
          REQUIRE(d_baseline.prec[k] == d_test.prec[k]);
          REQUIRE(d_baseline.cape[k] == Approx(d_test.cape[k]).margin(margin));
          REQUIRE(d_baseline.dcape[k] == Approx(d_test.dcape[k]).margin(margin));
          REQUIRE(d_baseline.dsubcld[k] == d_test.dsubcld[k]);
          REQUIRE(d_baseline.rliq[k] == d_test.rliq[k]);
          REQUIRE(d_baseline.gather_index[k] == d_test.gather_index[k]);
          REQUIRE(d_baseline.msemax_klev[k] == d_test.msemax_klev[k]);
          REQUIRE(d_baseline.jctop[k] == d_test.jctop[k]);
          REQUIRE(d_baseline.jcbot[k] == d_test.jcbot[k]);
          REQUIRE(d_baseline.jt[k] == d_test.jt[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.heat));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.qtnd));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.zdu));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.mflx_up));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.entr_up));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.detr_up));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.mflx_dn));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.entr_dn));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.p_del));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.ql));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.rprd));
        REQUIRE(d_baseline.total(d_baseline.heat) == d_test.total(d_test.dlf));
        for (Int k = 0; k < d_baseline.total(d_baseline.heat); ++k) {
          REQUIRE(d_baseline.heat[k] == d_test.heat[k]);
          REQUIRE(d_baseline.qtnd[k] == d_test.qtnd[k]);
          REQUIRE(d_baseline.zdu[k] == d_test.zdu[k]);
          REQUIRE(d_baseline.mflx_up[k] == Approx(d_test.mflx_up[k]).margin(margin));
          REQUIRE(d_baseline.entr_up[k] == Approx(d_test.entr_up[k]).margin(margin));
          REQUIRE(d_baseline.detr_up[k] == Approx(d_test.detr_up[k]).margin(margin));
          REQUIRE(d_baseline.mflx_dn[k] == Approx(d_test.mflx_dn[k]).margin(margin));
          REQUIRE(d_baseline.entr_dn[k] == d_test.entr_dn[k]);
          REQUIRE(d_baseline.p_del[k] == d_test.p_del[k]);
          REQUIRE(d_baseline.ql[k] == Approx(d_test.ql[k]).margin(margin));
          REQUIRE(d_baseline.rprd[k] == Approx(d_test.rprd[k]).margin(margin));
          REQUIRE(d_baseline.dlf[k] == d_test.dlf[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.mcon) == d_test.total(d_test.mcon));
        REQUIRE(d_baseline.total(d_baseline.mcon) == d_test.total(d_test.pflx));
        for (Int k = 0; k < d_baseline.total(d_baseline.mcon); ++k) {
          REQUIRE(d_baseline.mcon[k] == d_test.mcon[k]);
          REQUIRE(d_baseline.pflx[k] == Approx(d_test.pflx[k]).margin(margin));
        }
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        test_data[i].write(Base::m_ofile);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace zm
} // namespace scream

namespace {

TEST_CASE("zm_conv_main_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmConvMain;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
