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

    // Set up inputs.
    ZmConvMainData baseline_data[] = {
      //             pcols, ncol, pver, pverp, time_step, is_first_step, lengath (output, set to zero)
      ZmConvMainData(    4,    4,   72,    73,       2.0,         true,  0),
      ZmConvMainData(    4,    4,   72,    73,       3.0,         false , 0),
      ZmConvMainData(    4,    4,  128,   129,       4.0,         true, 0),
      ZmConvMainData(    4,    4,  128,   129,       5.0,         false , 0),
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
    std::vector<std::vector<bool>> actives;
    for (auto& d : test_data) {
      if (this->m_baseline_action == GENERATE) {
        zm_conv_main_f(d);
      }
      else {
        auto active = zm_conv_main(d);
        actives.push_back(active);
      }
    }

    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmConvMainData& d_baseline = baseline_data[i];
        ZmConvMainData& d_test = test_data[i];
        const Int ncol = d_test.ncol;
        const Int pver = d_test.pver;
        const Int pverp = d_test.pverp;
        REQUIRE(d_baseline.lengath == d_test.lengath);
        REQUIRE(d_baseline.total(d_baseline.prec) == ncol);
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
        Int inactive_cnt = 0;
        for (Int k = 0; k < d_baseline.total(d_baseline.prec); ++k) {
          const bool active_col = actives[i][k];
          if (!active_col) ++inactive_cnt;
          const Int fgindex = k-inactive_cnt;
          REQUIRE(d_baseline.prec[k] == Approx(d_test.prec[k]).margin(margin));
          REQUIRE(d_baseline.cape[k] == Approx(d_test.cape[k]).margin(margin));
          REQUIRE(d_baseline.dcape[k] == Approx(d_test.dcape[k]).margin(margin));
          REQUIRE(d_baseline.rliq[k] == Approx(d_test.rliq[k]).margin(margin));
          // Gathered 1-d variables
          if (active_col) {
            REQUIRE(d_baseline.dsubcld[fgindex] == d_test.dsubcld[k]);
            REQUIRE(d_baseline.msemax_klev[fgindex] == d_test.msemax_klev[k]);
            REQUIRE(d_baseline.jt[fgindex] == d_test.jt[k]);
            REQUIRE(d_baseline.jctop[k] == d_test.jctop[k]);
            REQUIRE(d_baseline.jcbot[k] == d_test.jcbot[k]);
          }
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
        REQUIRE(d_baseline.total(d_baseline.mcon) == d_test.total(d_test.mcon));
        REQUIRE(d_baseline.total(d_baseline.mcon) == d_test.total(d_test.pflx));
        inactive_cnt = 0;
        for (Int n = 0; n < ncol; ++n) {
          const bool active_col = actives[i][n];
          if (!active_col) ++inactive_cnt;
          for (Int k = 0; k < pverp; ++k) {
            const Int offset    = n*pver + k;
            const Int offsetp   = n*pverp + k;
            const Int fgoffset  = (n-inactive_cnt)*pver + k;
            const Int fgoffsetp = (n-inactive_cnt)*pverp + k;
            if (k < pver) {
              REQUIRE(d_baseline.mflx_up[offset] == Approx(d_test.mflx_up[offset]).margin(margin));
              REQUIRE(d_baseline.entr_up[offset] == Approx(d_test.entr_up[offset]).margin(margin));
              REQUIRE(d_baseline.detr_up[offset] == Approx(d_test.detr_up[offset]).margin(margin));
              REQUIRE(d_baseline.mflx_dn[offset] == Approx(d_test.mflx_dn[offset]).margin(margin));
              REQUIRE(d_baseline.entr_dn[offset] == d_test.entr_dn[offset]);
              // gathered variables: qtnd, rprd, zdu, mcon, heat, dlf, pflx, ql
              if (active_col) {
                REQUIRE(d_baseline.p_del[fgoffset] == d_test.p_del[offset]);
                REQUIRE(d_baseline.heat[offset] == d_test.heat[offset]);
                REQUIRE(d_baseline.qtnd[offset] == d_test.qtnd[offset]);
                REQUIRE(d_baseline.zdu[offset] == d_test.zdu[offset]);
                REQUIRE(d_baseline.ql[offset] == Approx(d_test.ql[offset]).margin(margin));
                REQUIRE(d_baseline.rprd[offset] == Approx(d_test.rprd[offset]).margin(margin));
                REQUIRE(d_baseline.dlf[offset] == d_test.dlf[offset]);
              }
            }
            // pverp variables (both gathered)
            if (active_col) {
              REQUIRE(d_baseline.mcon[offsetp] == d_test.mcon[offsetp]);
              REQUIRE(d_baseline.pflx[offsetp] == Approx(d_test.pflx[offsetp]).margin(margin));
            }
          }
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
