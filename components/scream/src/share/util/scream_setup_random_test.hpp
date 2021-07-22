#ifndef SCREAM_SETUP_RANDOM_TEST_HPP
#define SCREAM_SETUP_RANDOM_TEST_HPP

#include "catch2/catch.hpp"
#include "ekat/mpi/ekat_comm.hpp"

#include <random>
#include <iostream>

namespace scream {

/*
 * Create and return (via copy)
 */
template <typename Engine=std::mt19937_64>
Engine setup_random_test(const ekat::Comm* comm=nullptr)
{
  const auto& test_name = Catch::getResultCapture().getCurrentTestName();

  std::random_device rdev;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rdev() : catchRngSeed;

  if (comm == nullptr || comm->am_i_root()) {
    // Print seed to screen to trace tests that fail.
    std::cout << " For test " << test_name << ", random number generator seed: " << seed << "\n";
    if (catchRngSeed==0) {
      std::cout << "    Note: catch rng seed was 0 (default). We interpret that as a request to pick a random seed.\n"
        "    To reproduce a previous run, use --rng-seed N to provide the rng seed.\n\n";
    }
  }
  Engine engine(seed);
  return engine;
}

} // namespace scream

#endif // SCREAM_SETUP_RANDOM_TEST_HPP
