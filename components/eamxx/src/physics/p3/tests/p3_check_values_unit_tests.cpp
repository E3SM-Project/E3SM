#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

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
  auto engine = setup_random_test();

  CheckValuesData cvd_fortran[] = {
    //          kts_, kte_, timestepcount_, source_ind_, force_abort_
    CheckValuesData(1,  72,              2,         100,       false),
    CheckValuesData(1,  72,              3,         100,       false),
    CheckValuesData(1,  72,              4,         100,       false),
    CheckValuesData(1,  72,              5,         100,       false),
  };

  static constexpr Int num_runs = sizeof(cvd_fortran) / sizeof(CheckValuesData);

  for (auto& d : cvd_fortran) {
    d.randomize(engine, { {d.qv, {-4.056E-01, 1.153E+00}}, {d.temp, {1.000E+02, 5.000E+02}} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  CheckValuesData cvd_cxx[num_runs] = {
    CheckValuesData(cvd_fortran[0]),
    CheckValuesData(cvd_fortran[1]),
    CheckValuesData(cvd_fortran[2]),
    CheckValuesData(cvd_fortran[3]),
  };

  // Get data from fortran
  for (auto& d : cvd_fortran) {
    check_values(d);
  }

  // Get data from cxx
  for (auto& d : cvd_cxx) {
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
