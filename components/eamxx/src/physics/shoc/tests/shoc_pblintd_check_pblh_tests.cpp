#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPblintdCheckPblh : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
    static constexpr Int shcol = 5;
    static constexpr Int nlev = 4;
    static constexpr Int nlevi = nlev+1;

    // Tests for the subroutine pblintd_check_pblh

    // TEST ONE
    // Simple function and a simple test to be sure that "undefined"
    //   input values of PBLH are modified by this routine.

    // Define mid point height [m]
    static constexpr Real z[nlev] = {1500, 1000, 500, 20};
    // Define surface friction velocity [m/s]
    static constexpr Real ustar[shcol] = {0.1, 4.0, 0.9, 2.0, ustar_min};
    // Define check array
    static constexpr bool check[shcol] = {true, false, true, false, true};
    // Define PBL depth [m]
    // Make some values purposefully "undefined"
    static constexpr Real pblh[shcol] = {-20, 500, 1000, -20, 700};

    // Initialize data structure for bridging to F90
    PblintdCheckPblhData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(shcol > 0);
    REQUIRE(nlevi == nlev+1);

    // Fill in test data
    for(Int s = 0; s < shcol; ++s) {
      SDS.ustar[s] = ustar[s];
      SDS.check[s] = check[s];
      SDS.pblh[s] = pblh[s];
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.z[offset] = z[n];
      }
    }

    // check to make sure the input data makes sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.ustar[s] >= ustar_min);
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.z[offset] > 0);
      }
    }

    // Call the C++ implementation
    pblintd_check_pblh(SDS);

    // Check the result
    // Check that PBL height is greater than zero.  This is an
    //  important check to determine if the undefined values were modified.
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > 0);
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    PblintdCheckPblhData baseline_data[] = {
      PblintdCheckPblhData(36,  72, 73),
      PblintdCheckPblhData(72,  72, 73),
      PblintdCheckPblhData(128, 72, 73),
      PblintdCheckPblhData(256, 72, 73),
    };

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine, { {d.check, {1, 1}} });
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    PblintdCheckPblhData cxx_data[] = {
      PblintdCheckPblhData(baseline_data[0]),
      PblintdCheckPblhData(baseline_data[1]),
      PblintdCheckPblhData(baseline_data[2]),
      PblintdCheckPblhData(baseline_data[3]),
    };

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      pblintd_check_pblh(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(PblintdCheckPblhData);
      for (Int i = 0; i < num_runs; ++i) {
        PblintdCheckPblhData& d_baseline = baseline_data[i];
        PblintdCheckPblhData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.pblh); ++k) {
          REQUIRE(d_baseline.total(d_baseline.pblh) == d_cxx.total(d_cxx.pblh));
          REQUIRE(d_baseline.pblh[k] == d_cxx.pblh[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("pblintd_check_pblh_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCheckPblh;

  TestStruct().run_property();
}

TEST_CASE("pblintd_check_pblh_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCheckPblh;

  TestStruct().run_bfb();
}

} // empty namespace
