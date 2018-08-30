#include "catch2/catch.hpp"

#include "scream_util.hpp"

namespace {

TEST_CASE("precision", "util") {
  int nerr = 0;
  if(scream::util::is_single_precision<double>::value) nerr++;
  if(!scream::util::is_single_precision<float>::value) nerr++;
  REQUIRE(nerr == 0);
}

}  // namespace
