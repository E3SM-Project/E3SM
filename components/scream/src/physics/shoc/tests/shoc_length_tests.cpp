#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocLength {

  static void run_property()
  {
    static constexpr Real minlen = scream::shoc::Constants<Real>::minlen;
    static constexpr Real maxlen = scream::shoc::Constants<Real>::maxlen;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr Int nlevi    = nlev+1;

    // Tests for the upper level SHOC function:
    //   shoc_length

    // TEST ONE
    // Standard input and TKE test.
    // With reasonable inputs, verify outputs follow.
    // In addition, feed columns higher values of TKE and the length
    // scale for these columns should be gradually larger, given all other
    // inputs are equal.

    // PBL height [m]
    static constexpr Real pblh = 1000;
    // Define the host grid box size x-direction [m]
    static constexpr Real host_dx = 3000;
    // Defin the host grid box size y-direction [m]
    static constexpr Real host_dy = 5000;
    // Define the heights on the zt grid [m]
    static constexpr Real zi_grid[nlevi] = {9000, 5000, 1500, 900, 500, 0};
    // Virtual potential temperature on thermo grid [K]
    static constexpr Real thv[nlev] = {315, 310, 305, 300, 295};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec[nlev] = {0.02, 0.01, 0.04, 0.02, 0.05};
    // Turbulent kinetc energy [m2/s2]
    static constexpr Real tke[nlev] = {0.1, 0.15, 0.2, 0.25, 0.3};

    // compute geometric grid mesh
    const auto grid_mesh = sqrt(host_dx*host_dy);

    // set bound to check output
    static constexpr Real brunt_bound = 1;

    // Grid stuff to compute based on zi_grid
    Real zt_grid[nlev];
    Real dz_zt[nlev];
    Real dz_zi[nlevi];
    // Compute heights on midpoint grid
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      dz_zt[n] = zi_grid[n] - zi_grid[n+1];
      if (n == 0){
        dz_zi[n] = 0;
      }
      else{
        dz_zi[n] = zt_grid[n-1] - zt_grid[n];
      }
    }
    // set upper condition for dz_zi
    dz_zi[nlevi-1] = zt_grid[nlev-1];

    // Initialize data structure for bridging to F90
    ShocLengthData SDS(shcol, nlev, nlevi);

    // Load up input data
    for(Int s = 0; s < shcol; ++s) {
      SDS.host_dx[s] = host_dx;
      SDS.host_dy[s] = host_dy;
      SDS.pblh[s] = pblh;
      // Fill in test data on zt_grid.
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // for subsequent columns, increase TKE
        SDS.tke[offset] = (s+1)*tke[n];

        SDS.zt_grid[offset] = zt_grid[n];
        SDS.thv[offset] = thv[n];
        SDS.wthv_sec[offset] = wthv_sec[n];
        SDS.dz_zt[offset] = dz_zt[n];
      }

      // Fill in test data on zi_grid
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.zi_grid[offset] = zi_grid[n];
        SDS.dz_zi[offset] = dz_zi[n];
      }
    }

    // Check input data
    // for this test we need more than one column
    REQUIRE(SDS.shcol > 1);
    REQUIRE(SDS.nlevi == SDS.nlev+1);

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.zt_grid[offset] > 0);
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.thv[offset] > 0);
        REQUIRE(SDS.dz_zt[offset] > 0);

        // Make sure TKE is larger in next column over
        if (s < shcol-1){
          // get offset for "neighboring" column
          const auto offsets = n + (s+1) * nlev;
          REQUIRE(SDS.tke[offsets] > SDS.tke[offset]);
        }
      }

      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.zi_grid[offset] >= 0);
        REQUIRE(SDS.dz_zi[offset] >= 0);
      }
    }

    // Call the Fortran implementation
    shoc_length(SDS);

    // Verify output
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Require mixing length is greater than zero and is
        //  less than geometric grid mesh length + 1 m
        REQUIRE(SDS.shoc_mix[offset] >= minlen);
        REQUIRE(SDS.shoc_mix[offset] <= maxlen);
        REQUIRE(SDS.shoc_mix[offset] < 1.0+grid_mesh);

        // Be sure brunt vaisalla frequency is reasonable
        REQUIRE(std::abs(SDS.brunt[offset]) < 1);

        // Make sure length scale is larger when TKE is larger
        if (s < shcol-1){
          // get offset for "neighboring" column
          const auto offsets = n + (s+1) * nlev;
          if (SDS.tke[offsets] > SDS.tke[offset]){
            REQUIRE(SDS.shoc_mix[offsets] > SDS.shoc_mix[offset]);
          }
          else{
            REQUIRE(SDS.shoc_mix[offsets] < SDS.shoc_mix[offset]);
          }
        }
      }
    }

    // TEST TWO
    // Small grid mesh test.  Given a very small grid mesh, verify that
    //  the length scale is confined to this value.  Input from first
    //  test will be recycled into this one.

    // Define the host grid box size x-direction [m]
    static constexpr Real host_dx_small = 3;
    // Defin the host grid box size y-direction [m]
    static constexpr Real host_dy_small = 5;

    // compute geometric grid mesh
    const auto grid_mesh_small = sqrt(host_dx_small*host_dy_small);

    // Load new input
    for(Int s = 0; s < shcol; ++s) {
      SDS.host_dx[s] = host_dx_small;
      SDS.host_dy[s] = host_dy_small;
    }

    // call fortran implentation
    shoc_length(SDS);

    // Verify output
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Require mixing length is greater than zero and is
        //  less than geometric grid mesh length + 1 m
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.shoc_mix[offset] <= maxlen);
        REQUIRE(SDS.shoc_mix[offset] < 1.0+grid_mesh_small);
      }
    }

  }

  static void run_bfb()
  {
    ShocLengthData SDS_f90[] = {
      //        shcol, nlev, nlevi
      ShocLengthData(12, 71, 72),
      ShocLengthData(10, 12, 13),
      ShocLengthData(7,  16, 17),
      ShocLengthData(2, 7, 8),
    };

    static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(ShocLengthData);

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ShocLengthData SDS_cxx[] = {
      ShocLengthData(SDS_f90[0]),
      ShocLengthData(SDS_f90[1]),
      ShocLengthData(SDS_f90[2]),
      ShocLengthData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      shoc_length(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      shoc_length_f(d.shcol,d.nlev,d.nlevi,d.host_dx,d.host_dy,d.pblh,
                    d.tke,d.zt_grid,d.zi_grid,d.dz_zt,d.dz_zi,d.wthv_sec,d.thetal,
                    d.thv,d.brunt,d.shoc_mix);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      ShocLengthData& d_f90 = SDS_f90[i];
      ShocLengthData& d_cxx = SDS_cxx[i];
      for (Int k = 0; k < d_f90.total(d_f90.brunt); ++k) {
        REQUIRE(d_f90.brunt[k] == d_cxx.brunt[k]);
        REQUIRE(d_f90.shoc_mix[k] == d_cxx.shoc_mix[k]);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLength;

  TestStruct::run_property();
}

TEST_CASE("shoc_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLength;

  TestStruct::run_bfb();
}

} // namespace
