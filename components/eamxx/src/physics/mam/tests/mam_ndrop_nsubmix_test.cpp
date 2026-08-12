#include <catch2/catch.hpp>

#include <physics/mam/eamxx_mam_aci_process_interface.hpp>
#include <physics/mam/eamxx_mam_aci_functions.hpp>

TEST_CASE("ndrop nsubmix aggregation", "[mam4][ndrop]") {
  scream::Real latest = 100;
  scream::Real sum = 100;
  scream::Real maximum = 100;

  scream::detail::accumulate_ndrop_nsubmix(2, true, latest, sum, maximum);
  scream::detail::accumulate_ndrop_nsubmix(4, false, latest, sum, maximum);
  scream::detail::accumulate_ndrop_nsubmix(3, false, latest, sum, maximum);

  REQUIRE(latest == 3);
  REQUIRE(sum == 9);
  REQUIRE(maximum == 4);

  REQUIRE(scream::detail::can_accumulate_ndrop_nsubmix(
      3, false, scream::detail::max_exact_ndrop_nsubmix_aggregate - 3));
  REQUIRE_FALSE(scream::detail::can_accumulate_ndrop_nsubmix(
      4, false, scream::detail::max_exact_ndrop_nsubmix_aggregate - 3));
  REQUIRE_FALSE(scream::detail::can_accumulate_ndrop_nsubmix(0, true, 0));
  REQUIRE_FALSE(scream::detail::can_accumulate_ndrop_nsubmix(
      16777217, true, 0));

  sum = scream::detail::max_exact_ndrop_nsubmix_aggregate - 3;
  scream::detail::accumulate_ndrop_nsubmix(3, false, latest, sum, maximum);
  REQUIRE(latest == 3);
  REQUIRE(sum == scream::detail::max_exact_ndrop_nsubmix_aggregate);
  REQUIRE(maximum == 4);
}
