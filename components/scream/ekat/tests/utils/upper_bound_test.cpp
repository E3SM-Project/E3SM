#include <catch2/catch.hpp>

#include "ekat/util/scream_upper_bound.hpp"

#include <random>
#include <vector>
#include <algorithm>

namespace {

TEST_CASE("upper_bound", "soak") {
  std::default_random_engine generator;
  std::uniform_int_distribution<int> size_dist(100,1000);
  std::uniform_real_distribution<double> value_dist(0.0,1.0);

  for (int r = 0; r < 1000; ++r) {
    const int size = size_dist(generator);
    std::vector<double> v(size);
    for (int i = 0; i < size; ++i) {
      v[i] = value_dist(generator);
    }
    std::sort(v.begin(), v.end());
    double search_val = value_dist(generator);
    double v1 = *scream::util::upper_bound_impl(v.data(), v.data() + size, search_val);
    double v2 = *std::upper_bound(v.begin(), v.end(), search_val);
    REQUIRE(v1 == v2);
  }
}

} // empty namespace
