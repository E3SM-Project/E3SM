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
  CloudRainAccretionData crads_fortran[] = {
    // rho, inv_rho, qc_incld, nc_incld, qr_incld, qcacc, ncacc
    {1.0,  1.0, 72.0,   27.0,   72.0, 0.0, 0.0},
    {1.0,  1.0, 72.0,   72.0,   27.0, 0.0, 0.0},
    {1.0,  1.0, 72.0,   27.0,   72.0, 0.0, 0.0},
    {1.0,  1.0, 72.0,   72.0,   27.0, 0.0, 0.0},
    {1.0,  1.0, 72.0,   27.0,   27.0, 0.0, 0.0}
  };

  static constexpr Int num_runs = sizeof(crads_fortran) / sizeof(CloudRainAccretionData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  CloudRainAccretionData crads_cxx[num_runs] = {
    CloudRainAccretionData(crads_fortran[0]),
    CloudRainAccretionData(crads_fortran[1]),
    CloudRainAccretionData(crads_fortran[2]),
    CloudRainAccretionData(crads_fortran[3]),
    CloudRainAccretionData(crads_fortran[4]),
  };

  // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    cloud_rain_accretion(crads_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    CloudRainAccretionData& d = crads_cxx[i];
    cloud_rain_accretion_f(d.rho, d.inv_rho, d.qc_incld, d.nc_incld, d.qr_incld,
                           &d.qcacc, &d.ncacc);
  }

  // Check output.
  for (Int i = 0; i < num_runs; ++i) {
    REQUIRE(crads_fortran[i].qcacc == crads_cxx[i].qcacc);
    REQUIRE(crads_fortran[i].ncacc == crads_cxx[i].ncacc);
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
