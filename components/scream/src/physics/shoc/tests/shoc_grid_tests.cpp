#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_arch.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "physics/common/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

TEST_CASE("shoc_grid", "shoc") {
  constexpr Real gravit  = scream::physics::Constants<Real>::gravit;
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 128;
  constexpr auto nlevi   = nlev + 1;
  constexpr Real density = 1.0;
  constexpr Real dz      = 50.0;

  SHOCGridData SDS(shcol, nlev, nlevi);

  // Test that the inputs are resonable.
  REQUIRE(SDS.nlevi - SDS.nlev == 1);
  REQUIRE(SDS.shcol > 0);

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto nft    = SDS.nlev - 1 - n;
      const auto offset = n + s * SDS.nlev;

      SDS.zt_grid[offset] = nft * dz + dz / 2;
      SDS.pdel[offset]    = density * gravit * dz;
    }

    // Fill in test data on zi_grid.
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto nft      = SDS.nlevi - 1 - n;
      const auto offset   = n + s * SDS.nlevi;
      SDS.zi_grid[offset] = dz * nft;
    }
  }

  // Check that the inputs make sense

  // Check that zt decreases upward
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev - 1; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);
    }

    // Check that zi decreases upward
    for(Int n = 0; n < SDS.nlevi - 1; ++n) {
      const auto offset = n + s * SDS.nlevi;
      REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0.0);
    }
  }

  // Call the fortran implementation
  shoc_grid(nlev, SDS);

  // First check that dz is correct
  for(Int s = 0; s < shcol; ++s) {
    Real zt_sum = 0;
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.dz_zt[offset] > 0.0);
      REQUIRE(SDS.dz_zt[offset] == dz);
      zt_sum += SDS.dz_zt[offset];
    }
    // Check that the sum of dz_zt is equal to the largest zi
    REQUIRE(zt_sum == SDS.zi_grid[0]);
  }

  for(Int s = 0; s < shcol; ++s) {
    const auto s_offset = s * SDS.nlevi;
    REQUIRE(SDS.dz_zi[s_offset] == 0);
    REQUIRE(SDS.dz_zi[s_offset + SDS.nlevi - 1] == dz / 2.0);

    Real zi_sum = 0;
    for(Int n = 1; n < SDS.nlevi - 1; ++n) {
      const auto offset = n + s * SDS.nlevi;
      REQUIRE(SDS.dz_zi[offset] > 0.0);
      REQUIRE(SDS.dz_zi[offset] == dz);
      zi_sum += SDS.dz_zi[offset];
    }
    // Check that the sum of dz_zi is equal to the largest zt
    zi_sum += SDS.dz_zi[SDS.nlevi - 1];
    REQUIRE(zi_sum == SDS.zt_grid[0]);
  }

  // Now check density
  for(Int s = 0; s < shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      // check that the density is consistent with the hydrostatic approximation
      REQUIRE(SDS.rho_zt[offset] == density);

      // check that the density has physically realistic values
      REQUIRE(SDS.rho_zt[offset] <= 2.0);
      REQUIRE(SDS.rho_zt[offset] > 0.0);
    }
  }
}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
