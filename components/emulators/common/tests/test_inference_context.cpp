// Catch2 v2 single header
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "inference_context.hpp"
#include "inference_error.hpp"

#include <string>

namespace emulator {
namespace inference {
namespace test {

TEST_CASE("A default context describes a serial run", "[context]") {
  InferenceContext context;
  REQUIRE(context.rank == 0);
  REQUIRE(context.size == 1);
  REQUIRE(context.is_root());
  REQUIRE(context.num_local_cols() == 0);
}

TEST_CASE("The context carries the coupler's decomposition", "[context]") {
  const int gids[3] = {1, 5, 9};
  const double lat[3] = {-45.0, 0.0, 45.0};
  const double lon[3] = {0.0, 120.0, 240.0};

  InferenceContext context;
  context.set_grid(8, 6, 48, gids, lat, lon, 3);

  REQUIRE(context.num_local_cols() == 3);
  REQUIRE(context.col_gids[2] == 9);
  REQUIRE(context.lat[0] == -45.0);
  REQUIRE(context.nx == 8);
  REQUIRE(context.to_string().find("3 of 48 columns") != std::string::npos);
}

TEST_CASE("make_context is usable without a launcher", "[context]") {
  // Built with MPI but never launched under it (a unit test, a tool), the
  // context must still describe a valid serial run rather than abort.
  const auto context = make_context(0);
  REQUIRE(context.size >= 1);
  REQUIRE(context.rank >= 0);
  REQUIRE(context.rank < context.size);
}

} // namespace test
} // namespace inference
} // namespace emulator
