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
struct UnitWrap::UnitTest<D>::TestPblintdSurfTemp : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
    static constexpr Int shcol = 4;
    static constexpr Int nlev = 4;
    static constexpr Int nlevi = nlev+1;

    // Tests for the subroutine pblintd_surf_temp

    // Define mid point height [m]
    static constexpr Real z[nlev] = {1500, 1000, 500, 20};
    // Define virtual potential temperature [K]
    static constexpr Real thv[nlev]= {305, 302, 302, 300};
    // Define obklen [m]
    static constexpr Real obklen[shcol] = {-10, 20, -100, 150};
    // Define surf friction velocity [m4/s4]
    static constexpr Real ustar[shcol] = {ustar_min, 1, 5, 10};
    // Surface kinematic buoyancy flux
    static constexpr Real kbfs[shcol] = {0.03, -0.03, 0.1, -0.1};
    // Input check parameter
    static constexpr bool check[shcol] = {true, true, false, false};
    // Input bulk richardson number
    static constexpr Real rino = 0.1;

    // Initialize data structure for bridging to F90
    PblintdSurfTempData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(shcol > 0);
    REQUIRE(nlevi == nlev+1);

    // Fill in test data on the zt_grid
    for(Int s = 0; s < shcol; ++s) {
      SDS.ustar[s] = ustar[s];
      SDS.check[s] = check[s];
      SDS.kbfs[s] = kbfs[s];
      SDS.obklen[s] = obklen[s];
      SDS.pblh[s] = 0; // Initialize to zero
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.z[offset] = z[n];
        SDS.thv[offset] = thv[n];
        SDS.rino[offset] = rino;
      }
    }

    // check to make sure the input data makes sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.ustar[s] >= ustar_min);
      REQUIRE(std::abs(SDS.kbfs[s]) < 1);
      // Make sure the sign of the Monin Obukov length is
      //  consistent with the sign of kinematic buoyancy flux
      if (SDS.kbfs[s] < 0){
        REQUIRE(SDS.obklen[s] > 0);
      }
      else{
        REQUIRE(SDS.obklen[s] < 0);
      }
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.z[offset] > 0);
        REQUIRE(SDS.thv[offset] < 350);
        REQUIRE(SDS.thv[offset] > 290);
      }
    }

    // Call the C++ implementation
    pblintd_surf_temp(SDS);

    // Check the result
    for(Int s = 0; s < shcol; ++s) {
      // If OUTPUT check was performed, verify that surface diagnosed
      //  temperature is appropriate
      if (SDS.check[s] == true){
        REQUIRE(SDS.tlv[s] < 350);
        REQUIRE(SDS.tlv[s] > 290);

       // determine lowest level index for this column
       const auto low_lev = (nlev-1) + s * nlev;
       REQUIRE(SDS.rino[low_lev] == 0);
      }

      // If INPUT check was performed, verify that PBLH depth was
      //  set to the expected value, which should be the highest level
      //  of the data provided
      if (check[s] == true){
        REQUIRE(SDS.pblh[s] == SDS.z[0]);
      }
    }

  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    PblintdSurfTempData baseline_data[] = {
      PblintdSurfTempData(6, 7, 8),
      PblintdSurfTempData(64, 72, 73),
      PblintdSurfTempData(128, 72, 73),
      PblintdSurfTempData(256, 72, 73),
    };

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine, { {d.obklen, {100., 200.}} });
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    PblintdSurfTempData cxx_data[] = {
      PblintdSurfTempData(baseline_data[0]),
      PblintdSurfTempData(baseline_data[1]),
      PblintdSurfTempData(baseline_data[2]),
      PblintdSurfTempData(baseline_data[3]),
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
      pblintd_surf_temp(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(PblintdSurfTempData);
      for (Int i = 0; i < num_runs; ++i) {
        PblintdSurfTempData& d_baseline = baseline_data[i];
        PblintdSurfTempData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tlv); ++k) {
          REQUIRE(d_baseline.total(d_baseline.tlv) == d_cxx.total(d_cxx.tlv));
          REQUIRE(d_baseline.tlv[k] == d_cxx.tlv[k]);
          REQUIRE(d_baseline.total(d_baseline.tlv) == d_cxx.total(d_cxx.pblh));
          REQUIRE(d_baseline.pblh[k] == d_cxx.pblh[k]);
          REQUIRE(d_baseline.total(d_baseline.tlv) == d_cxx.total(d_cxx.check));
          REQUIRE(d_baseline.check[k] == d_cxx.check[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.rino); ++k) {
          REQUIRE(d_baseline.total(d_baseline.rino) == d_cxx.total(d_cxx.rino));
          REQUIRE(d_baseline.rino[k] == d_cxx.rino[k]);
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

TEST_CASE("pblintd_surf_temp_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdSurfTemp;

  TestStruct().run_property();
}

TEST_CASE("pblintd_surf_temp_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdSurfTemp;

  TestStruct().run_bfb();
}

} // empty namespace
