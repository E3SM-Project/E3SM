#include <catch2/catch.hpp>

#include "moab_adapter.hpp"

#include <mpi.h>

#include <span>
#include <string>
#include <vector>

using namespace e3sm::emulator::coupler;

TEST_CASE("MoabAdapter stores scalar and vector fields", "[moab][mpi]") {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  INFO("MPI rank: " << rank);

  constexpr std::size_t nlocal = 2;
  constexpr std::size_t nlevels = 3;
  constexpr double sentinel = -999.0;

  // Unique IDs across the component communicator.
  std::vector<int> gids{
      rank * 2 + 1,
      rank * 2 + 2,
  };

  std::vector<double> latitudes{10.0, 20.0};
  std::vector<double> longitudes{
      -100.0 + rank,
      -99.5 + rank,
  };

  MoabPartition partition{};
  partition.comm = MPI_COMM_WORLD;
  partition.global_ids = std::span<int>{gids};
  partition.latitudes = std::span<double>{latitudes};
  partition.longitudes = std::span<double>{longitudes};

  // Same logical component name and ID on every participating rank.
  MoabAdapter adapter{"adapter-test", 1, partition};

  adapter.define_tag(MoabTagDescriptor{
      .name = "temperature",
      .role = FieldRole::Export,
      .values_per_vertex = 1,
  });

  adapter.define_tag(MoabTagDescriptor{
      .name = "soil_moisture",
      .role = FieldRole::Export,
      .values_per_vertex = nlevels,
  });

  std::vector<double> temperature(nlocal);
  std::vector<double> moisture(nlocal * nlevels);

  for (std::size_t column = 0; column < nlocal; ++column) {
    temperature[column] = 280.0 + gids[column];

    for (std::size_t level = 0; level < nlevels; ++level) {
      moisture[column * nlevels + level] = 0.01 * gids[column] + 0.001 * level;
    }
  }

  adapter.write_tag("temperature", std::span<double>{temperature});
  adapter.write_tag("soil_moisture", std::span<double>{moisture});

  std::vector<double> recovered_temperature(nlocal, sentinel);
  std::vector<double> recovered_moisture(nlocal * nlevels, sentinel);

  adapter.read_tag("temperature", std::span<double>{recovered_temperature});
  adapter.read_tag("soil_moisture", std::span<double>{recovered_moisture});

  // CHECK continues execution on all ranks if a comparison fails.
  REQUIRE(recovered_temperature == temperature);
  REQUIRE(recovered_moisture == moisture);

  // Changing caller-owned data must not change the stored tag values.
  const std::vector<double> original_temperature = temperature;
  temperature[0] += 10.0;

  adapter.read_tag("temperature", recovered_temperature);
  REQUIRE(recovered_temperature == original_temperature);

  // Explicitly writing again updates MOAB's stored values.
  adapter.write_tag("temperature", temperature);
  adapter.read_tag("temperature", recovered_temperature);

  REQUIRE(recovered_temperature == temperature);

  // These checks must fail locally, before entering any MOAB operation.
  std::vector<double> wrong_size(nlocal - 1, sentinel);

  CHECK_THROWS_AS(adapter.read_tag("temperature", wrong_size),
                  std::runtime_error);

  CHECK_THROWS_AS(adapter.write_tag("temperature", wrong_size),
                  std::runtime_error);

  CHECK_THROWS_AS(adapter.read_tag("missing_tag", recovered_temperature),
                  std::runtime_error);

  CHECK_THROWS_AS(adapter.write_tag("missing_tag", temperature),
                  std::runtime_error);

  // adapter is destroyed before the runner finalizes MOAB.
}
