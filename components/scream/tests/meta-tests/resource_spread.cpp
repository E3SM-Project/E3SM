#include "ekat/kokkos/ekat_kokkos_types.hpp"

#include <catch2/catch.hpp>

#include <iostream>
#include <chrono>
#include <thread>

namespace scream {

TEST_CASE("rank_and_thread_spread", "[fake_infra_test]")
{
  // Nothing needs to be done here except sleeping to give a chance
  // for tests to run concurrently.
  const auto seconds_to_sleep = 5;
  std::this_thread::sleep_for(std::chrono::seconds(seconds_to_sleep));
}

} // empty namespace
