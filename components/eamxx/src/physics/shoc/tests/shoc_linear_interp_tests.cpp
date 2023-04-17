#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocLinearInt {

  static void run_property_fixed()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int km1     = 5;
    static constexpr auto km2   = km1 + 1;

    // TEST ONE
    // Test interpolation going from midpoint grid to interface grid.  Note that
    //  in this case the nlev grid is denoted by variable km1 and the nlevi
    //  grid is represented by variable km2.  This is because the linear
    //  interp routine can interpolate from midpoint to interface grid or
    //  vice versa so notation should be flexible.

    // Note that in this test we are testing TWO profils.  The first column
    //  will be loaded up with potential temperature values while the second
    //  column will be loaded up with meridional wind values.

    // Define the interface height grid [m]
    static constexpr Real zi_grid[km2] = {12500., 7500., 3000., 750., 250.0, 0.};
    // Define the liquid water potential temperature [K]
    //  on the midpoint grid
    static constexpr Real thetal_zt[km1] = {320.0, 310.0, 300.0, 300.0, 306.0};
    // Define meridional wind on midpoint grid
    static constexpr Real u_wind_zt[km1] = {-2.0, -0.5, 0.0, 5.0, 2.0};
    // Define minimum threshold for this experiment, since we are
    //  dealing with negative winds in a column then set to a large
    //  negative value since we do not want to clip.
    static constexpr Real minthresh = -99999.0;

    // Initialize data structure for bridging to F90
    LinearInterpData SDS(shcol, km1, km2, minthresh);

    // For this test we need exactly two columns
    REQUIRE( (SDS.ncol == shcol && SDS.km1 == km1 && SDS.km2 == km2 && SDS.minthresh == minthresh) );
    REQUIRE(shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < km1; ++n) {
        const auto offset = n + s * km1;

        // For zt grid heights compute as midpoint
        //  between interface heights
        SDS.x1[offset] = 0.5*(zi_grid[n]+zi_grid[n+1]);
        if (s == 0){
          SDS.y1[offset] = thetal_zt[n];
        }
        else{
          SDS.y1[offset] = u_wind_zt[n];
        }
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < km2; ++n) {
        const auto offset   = n + s * km2;
        SDS.x2[offset] = zi_grid[n];
      }
    }

    // Check that the inputs make sense

    // Check that zt decreases upward
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < km1 - 1; ++n) {
        const auto offset = n + s * km1;
        REQUIRE(SDS.x1[offset + 1] - SDS.x1[offset] < 0.0);
      }

      // Check that zi decreases upward
      for(Int n = 0; n < km2 - 1; ++n) {
        const auto offset = n + s * km2;
        REQUIRE(SDS.x2[offset + 1] - SDS.x2[offset] < 0.0);
      }
    }

    // Call the fortran implementation
    linear_interp(SDS);

    // First check that all output temperatures are greater than zero

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < km1; ++n) {
        const auto offset = n + s * km1;
        const auto offseti = n + s * km2;

        // check boundary points
        // First upper boundary
        if (n == 0){
          const auto uppergradient = SDS.y1[offset] - SDS.y1[offset+1];
          if (uppergradient > 0.0){
            // if upper gradient is positive then make sure that
            //  the temperature of the highest nlevi layer is greater
            //  than the tempature of the highest nlev layer
            REQUIRE(SDS.y2[offseti] > SDS.y1[offset]);
          }
          if (uppergradient < 0.0){
            REQUIRE(SDS.y2[offseti] < SDS.y1[offset]);
          }
          if (uppergradient == 0.0){
            REQUIRE(SDS.y2[offseti] == SDS.y2[offseti]);
          }
        }

        // Now the lower boundary
        else if (n == km1-1){
          const auto lowergradient = SDS.y1[offset-1] - SDS.y1[offset];
          if (lowergradient > 0.0){
            // if lower gradient is positive then make sure that
            //  the temperature of the lowest nlevi layer is lower
            //  than the tempature of the lowest nlev layer
            REQUIRE(SDS.y2[offseti+1] < SDS.y1[offset]);
          }
          if (lowergradient < 0.0){
            REQUIRE(SDS.y2[offseti+1] > SDS.y1[offset]);
          }
          if (lowergradient == 0.0){
            REQUIRE(SDS.y2[offseti+1] == SDS.y2[offset]);
          }
        }

        // Now make sure all points are bounded as expected
        else{
          const auto gradient = SDS.y1[offset-1] - SDS.y1[offset];
          if (gradient == 0.0){
            REQUIRE(SDS.y2[offseti] == SDS.y1[offset]);
          }
          else if (gradient > 0.0){
            REQUIRE(SDS.y2[offseti] < SDS.y1[offset-1]);
            REQUIRE(SDS.y2[offseti] > SDS.y1[offset]);
          }
          else {
            REQUIRE(SDS.y2[offseti] > SDS.y1[offset-1]);
            REQUIRE(SDS.y2[offseti] < SDS.y1[offset]);
          }
        }

      }
    }

//  TEST TWO
//  Now we test going from interface grid to mid point grid

    // Initialize data structure for bridging to F90
    // NOTE that km2 and km1 grid must be switched here.
    // Must initialize a new data structure since km1 and km2 are swapped.
    LinearInterpData SDS2(shcol, km2, km1, minthresh);

    // Fill in test data on zi_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < km2; ++n) {
        const auto offset = n + s * km2;

        // Load up stuff on interface grid.  Here we
        //  are going to use information from the last test
        //  to initialize our grid
        SDS2.x1[offset] = SDS.x2[offset];
        SDS2.y1[offset] = SDS.y2[offset];
      }

      // Load up stuff on midpoint grid, use output from
      //  the last test here.
      for(Int n = 0; n < km1; ++n) {
        const auto offset   = n + s * km1;
        SDS2.x2[offset] = SDS.x1[offset];
      }
    }

    linear_interp(SDS2);

    // Check the result, make sure output is bounded correctly

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < km1; ++n) {
        const auto offset = n + s * km1;
        const auto offseti = n + s * km2;

        // compute gradient
        const auto gradient = SDS2.y1[offseti] - SDS2.y1[offseti+1];

        if (gradient == 0.0){
          REQUIRE(SDS2.y2[offset] == SDS2.y1[offseti]);
        }
        else if (gradient > 0.0){
          REQUIRE(SDS2.y2[offset] > SDS2.y1[offseti+1]);
          REQUIRE(SDS2.y2[offset] < SDS2.y1[offseti]);
        }
        else {
          REQUIRE(SDS2.y2[offset] < SDS2.y1[offseti+1]);
          REQUIRE(SDS2.y2[offset] > SDS2.y1[offseti]);
        }
      }
    }
  }

  static void run_property_random(bool km1_bigger)
  {
    std::default_random_engine generator;
    std::pair<Int, Int> km1_range = {13, 25};
    std::pair<Int, Int> shcol_range = {5, 10};
    std::pair<Real, Real> y1_range = {0, 1};
    std::pair<Real, Real> x_range = {0, 1};

    static constexpr Real minthresh = -99999.0;

    std::uniform_int_distribution<Int> km1_dist(km1_range.first, km1_range.second);
    std::uniform_int_distribution<Int> shcol_dist(shcol_range.first, shcol_range.second);
    std::uniform_real_distribution<Real> y1_dist(y1_range.first, y1_range.second);
    std::uniform_real_distribution<Real> x_dist(x_range.first, x_range.second);

    const Int shcol = shcol_dist(generator);
    const Int km1 = km1_dist(generator);
    const Int km2 = km1 + (km1_bigger ? -1 : 1);

    LinearInterpData d(shcol, km1, km2, minthresh);

    for (Int s = 0; s < shcol; ++s) {

      if (km1_bigger) {
        for (Int k1 = 0; k1 < km1; ++k1) {
          const auto offset = k1 + s*km1;
          d.x1[offset] = x_dist(generator);
          d.y1[offset] = y1_dist(generator);
        }
        std::sort(d.x1 + (s*km1), d.x1 + ((s+1)*km1));

        for (Int k2 = 0; k2 < km2; ++k2) {
          const auto offset  = k2 + s*km2;
          const auto offset1 = k2 + s*km1;
          REQUIRE(d.x1[offset1] < d.x1[offset1+1]);
          std::uniform_real_distribution<Real> x2_dist(d.x1[offset1], d.x1[offset1+1]);
          d.x2[offset] = x2_dist(generator);
        }
      }
      else {
        for (Int k2 = 0; k2 < km2; ++k2) {
          const auto offset = k2 + s*km2;
          d.x2[offset] = x_dist(generator);
        }
        std::sort(d.x2 + (s*km2), d.x2 + ((s+1)*km2));

        for (Int k1 = 0; k1 < km1; ++k1) {
          const auto offset  = k1 + s*km1;
          const auto offset2 = k1 + s*km2;
          REQUIRE(d.x2[offset2] < d.x2[offset2+1]);
          std::uniform_real_distribution<Real> x1_dist(d.x2[offset2], d.x2[offset2+1]);
          d.x1[offset] = x1_dist(generator);
          d.y1[offset] = y1_dist(generator);
        }
      }
    }

    linear_interp(d);

    // The combination of single-precision and randomness generating points
    // close together can result in larger error margins.
    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    for (Int s = 0; s < shcol; ++s) {
      if (km1_bigger) {
        for (Int k2 = 0; k2 < km2; ++k2) {
          const auto offset  = k2 + s*km2;
          const auto offset1 = k2 + s*km1;
          const auto x1_delta = d.x1[offset1+1] - d.x1[offset1];
          const auto y1_delta = d.y1[offset1+1] - d.y1[offset1];
          const auto slope1   = y1_delta / x1_delta;
          const auto x2_delta = d.x2[offset] - d.x1[offset1];
          const auto y2_delta = d.y2[offset] - d.y1[offset1];
          const auto slope2   = y2_delta / x2_delta;
          REQUIRE(slope1 == Approx(slope2).margin(margin));
        }
      }
      else {
        for (Int k1 = 0; k1 < km1; ++k1) {
          if (k1 == 0 || k1 == km1-1) {
            const auto offset  = k1 + s*km1;
            const auto offset2 = k1 + s*km2;
            const auto x2_delta = d.x2[offset2+1] - d.x2[offset2];
            const auto y2_delta = d.y2[offset2+1] - d.y2[offset2];
            const auto slope2   = y2_delta / x2_delta;
            const auto x1_delta = d.x1[offset] - d.x2[offset2];
            const auto y1_delta = d.y1[offset] - d.y2[offset2];
            const auto slope1   = y1_delta / x1_delta;
            REQUIRE(slope1 == Approx(slope2).margin(margin));
          }
          else {
            const auto offset  = k1 + s*km1;
            const auto offset2 = k1 + s*km2;
            const auto x1_delta = d.x1[offset] - d.x1[offset-1];
            const auto y1_delta = d.y1[offset] - d.y1[offset-1];
            const auto slope1   = y1_delta / x1_delta;
            const auto x2_delta = d.x2[offset2] - d.x1[offset-1];
            const auto y2_delta = d.y2[offset2] - d.y1[offset-1];
            const auto slope2   = y2_delta / x2_delta;
            REQUIRE(slope1 == Approx(slope2).margin(margin));
          }
        }
      }
    }
  }

  static void run_property()
  {
    run_property_fixed();
    run_property_random(true);
    run_property_random(false);
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    LinearInterpData f90_data[] = {
      //                   shcol, nlev(km1), nlevi(km2), minthresh
      LinearInterpData(10, 72, 71, 1e-15),
      LinearInterpData(10, 71, 72, 1e-15),
      LinearInterpData(1, 15, 16, 1e-15),
      LinearInterpData(1, 16, 15, 1e-15),
      LinearInterpData(1, 5, 6, 1e-15),
      LinearInterpData(1, 6, 5, 1e-15),
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    LinearInterpData cxx_data[] = {
      LinearInterpData(f90_data[0]),
      LinearInterpData(f90_data[1]),
      LinearInterpData(f90_data[2]),
      LinearInterpData(f90_data[3]),
      LinearInterpData(f90_data[4]),
      LinearInterpData(f90_data[5]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      linear_interp(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      linear_interp_f(d.x1, d.x2, d.y1, d.y2, d.km1, d.km2, d.ncol, d.minthresh);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(LinearInterpData);
      for (Int i = 0; i < num_runs; ++i) {
        LinearInterpData& d_f90 = f90_data[i];
        LinearInterpData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.y2); ++k) {
          REQUIRE(d_f90.y2[k] == d_cxx.y2[k]);
        }
      }
    }
  } // run_bfb

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_linear_interp_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLinearInt;

  TestStruct::run_property();
}

TEST_CASE("shoc_linear_interp_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLinearInt;

  TestStruct::run_bfb();
}

} // namespace
