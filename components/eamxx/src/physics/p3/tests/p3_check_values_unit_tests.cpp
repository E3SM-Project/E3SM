#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_test_data.hpp"
#include "p3_data.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestCheckValues : public UnitWrap::UnitTest<D>::Base {

void run_check_values_bfb()
{
  // This is not really a bfb test since no results are being checked.
  auto engine = setup_random_test();

  CheckValuesData cvd_cxx[] = {
    //          kts_, kte_, timestepcount_, source_ind_, force_abort_
    CheckValuesData(1,  72,              2,         100,       false),
    CheckValuesData(1,  72,              3,         100,       false),
    CheckValuesData(1,  72,              4,         100,       false),
    CheckValuesData(1,  72,              5,         100,       false),
  };

  for (auto& d : cvd_cxx) {
    d.randomize(engine, { {d.qv, {-4.056E-01, 1.153E+00}}, {d.temp, {1.000E+02, 5.000E+02}} });
  }

  // Get data from cxx
  for (auto& d : cvd_cxx) {
    check_values_host(d.qv, d.temp, d.kts, d.kte, d.timestepcount, d.force_abort, d.source_ind, d.col_loc);
  }
}

void run_check_values_phys()
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
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCheckValues;

  T t;
  t.run_check_values_phys();
  t.run_check_values_bfb();
}

} // namespace
