#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/p3/p3_f90.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestCheckValues {

static void run_check_values_bfb()
{
  const std::array< std::pair<Real, Real>, CheckValuesData::NUM_ARRAYS > ranges = {
    std::make_pair(-4.056E-01, 1.153E+00), // qv
    std::make_pair( 1.000E+02, 5.000E+02), // temp
  };

  CheckValuesData cvd_fortran[] = {
    // kts_, kte_, timestepcount_, source_ind_, force_abort_, ranges
    CheckValuesData(1,  72,   2,   100,   false,  ranges),
    CheckValuesData(1,  72,   3,   100,   false,  ranges),
    CheckValuesData(1,  72,   4,   100,   false,  ranges),
    CheckValuesData(1,  72,   5,   100,   false,  ranges),
  };

  static constexpr Int num_runs = sizeof(cvd_fortran) / sizeof(CheckValuesData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  CheckValuesData cvd_cxx[num_runs] = {
    CheckValuesData(cvd_fortran[0]),
    CheckValuesData(cvd_fortran[1]),
    CheckValuesData(cvd_fortran[2]),
    CheckValuesData(cvd_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    check_values(cvd_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    CheckValuesData& d = cvd_cxx[i];
    check_values_f(d.qv, d.temp, d.kts, d.kte, d.timestepcount, d.force_abort, d.source_ind, d.col_loc);
  }
}

static void run_check_values_phys()
{
    // TODO
}

};

}
}
}

namespace {

TEST_CASE("p3_check_values", "[p3_functions]")
{
  using TRS = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCheckValues;

  scream::p3::p3_init(); // need fortran table data

  TRS::run_check_values_phys();
  TRS::run_check_values_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
