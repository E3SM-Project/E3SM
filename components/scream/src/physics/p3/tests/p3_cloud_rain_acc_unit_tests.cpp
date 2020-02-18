#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestCloudRainAccretion {

static void run_phys()
{
  // TODO
}

static void run_bfb()
{
  // This is the threshold for whether the qc and qr cloud mixing ratios are
  // large enough to affect the warm-phase process rates qcacc and ncacc.
  constexpr Scalar qsmall = C::QSMALL;

  constexpr Scalar rho = 1.0;
  constexpr Scalar inv_rho = 1.0 / rho;
  constexpr Scalar qc_incld_small = 0.9 * qsmall;
  constexpr Scalar qr_incld_small = 0.9 * qsmall;
  constexpr Scalar qc_incld_not_small = 2.0 * qsmall;
  constexpr Scalar qr_incld_not_small = 2.0 * qsmall;
  constexpr Scalar nc_incld = 1.0;

  CloudRainAccretionData fortran_data[] = {
      // rho, inv_rho, qc_incld, nc_incld, qr_incld, qcacc, ncacc
      {rho, inv_rho, qc_incld_small, nc_incld, qr_incld_small, 0.0, 0.0},
      {rho, inv_rho, qc_incld_small, nc_incld, qr_incld_not_small, 0.0, 0.0},
      {rho, inv_rho, qc_incld_not_small, nc_incld, qr_incld_small, 0.0, 0.0},
      {rho, inv_rho, qc_incld_not_small, nc_incld, qr_incld_not_small, 0.0, 0.0},
  };

  static constexpr Int num_runs = sizeof(fortran_data) / sizeof(CloudRainAccretionData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  CloudRainAccretionData cxx_data[] = {
      fortran_data[0],
      fortran_data[1],
      fortran_data[2],
      fortran_data[3]
  };

  // Run the Fortran subroutine.
  for (Int i = 0; i < num_runs; ++i) {
    cloud_rain_accretion(fortran_data[i]);
  }

  // Run the C++ subroutine.
  for (Int i = 0; i < num_runs; ++i) {
    CloudRainAccretionData& d = cxx_data[i];
    cloud_rain_accretion_f(d.rho, d.inv_rho, d.qc_incld, d.nc_incld, d.qr_incld,
                           &d.qcacc, &d.ncacc);
  }

  // Check output.
  for (Int i = 0; i < num_runs; ++i) {
    REQUIRE(fortran_data[i].qcacc == cxx_data[i].qcacc);
    REQUIRE(fortran_data[i].ncacc == cxx_data[i].ncacc);
  }

}

};

}
}
}

namespace {

TEST_CASE("p3_cloud_rain_accretion", "[p3_functions]")
{
  using TCRA = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCloudRainAccretion;

  TCRA::run_phys();
  TCRA::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
