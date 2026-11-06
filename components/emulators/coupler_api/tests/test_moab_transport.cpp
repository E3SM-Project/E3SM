// test_moab_transport.cpp
#include <catch2/catch.hpp>

#include <mpi.h>

#include "moab_adapter.hpp"
#include "moab_transport.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <string>
#include <vector>

using namespace e3sm::coupler;
struct CommOwner {
  MPI_Comm value;
  ~CommOwner() {
    if (value != MPI_COMM_NULL) {
      MPI_Comm_free(&value);
    }
  }
};

TEST_CASE("Transport redistributes fields between component rank groups",
          "[moab][transport][mpi]") {
  int world_rank = 0;
  int world_size = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  INFO("World rank: " << world_rank);

  REQUIRE(world_size == 6);

  constexpr int global_columns = 12;
  constexpr std::size_t nlevels = 3;
  constexpr double sentinel = -999.0;

  const bool is_forcing = world_rank < 2;

  MPI_Comm component_comm = MPI_COMM_NULL;
  MPI_Comm_split(MPI_COMM_WORLD, is_forcing ? 0 : 1, world_rank,
                 &component_comm);

  // Cleanup helper!!
  CommOwner component_comm_owner{component_comm};

  int component_rank = 0;
  int component_size = 0;
  MPI_Comm_rank(component_comm, &component_rank);
  MPI_Comm_size(component_comm, &component_size);

  INFO("Component rank: " << component_rank);
  INFO("Component: " << (is_forcing ? "forcing" : "land"));

  const std::size_t nlocal =
      static_cast<std::size_t>(global_columns / component_size);

  std::vector<int> gids(nlocal);
  std::vector<double> latitudes(nlocal);
  std::vector<double> longitudes(nlocal);

  for (std::size_t local = 0; local < nlocal; ++local) {
    const int gid =
        component_rank * static_cast<int>(nlocal) + static_cast<int>(local) + 1;

    gids[local] = gid;
    latitudes[local] = 10.0 + gid;
    longitudes[local] = -100.0 + gid;
  }

  MoabPartition partition{};
  partition.comm = component_comm;
  partition.global_ids = std::span<int>{gids};
  partition.latitudes = std::span<double>{latitudes};
  partition.longitudes = std::span<double>{longitudes};

  // Declared before Transport, so the adapter outlives it.
  MoabAdapter adapter{
      is_forcing ? "forcing" : "land",
      is_forcing ? 1 : 2,
      partition,
  };

  const FieldRole role = is_forcing ? FieldRole::Export : FieldRole::Import;

  adapter.define_tag(MoabTagDescriptor{
      .name = "temperature",
      .role = role,
      .values_per_vertex = 1,
  });

  adapter.define_tag(MoabTagDescriptor{
      .name = "soil_moisture",
      .role = role,
      .values_per_vertex = nlevels,
  });

  std::array<MoabAdapter*, 1> local_adapters{&adapter};

  const std::array<TransportRoute, 2> routes{
      TransportRoute{
          .field = "temperature",
          .source_component = "forcing",
          .destination_component = "land",
      },
      TransportRoute{
          .field = "soil_moisture",
          .source_component = "forcing",
          .destination_component = "land",
      },
  };

  Transport transport{MPI_COMM_WORLD};
  transport.initialize(local_adapters, routes);

  // Integer-valued doubles make exact comparisons appropriate.
  const auto expected_temperature = [](int gid) {
    return 280.0 + static_cast<double>(gid);
  };

  const auto expected_moisture = [](int gid, std::size_t level) {
    return 100.0 * static_cast<double>(gid) + static_cast<double>(level);
  };

  std::vector<double> temperature(nlocal, sentinel);
  std::vector<double> moisture(nlocal * nlevels, sentinel);

  if (is_forcing) {
    for (std::size_t local = 0; local < nlocal; ++local) {
      temperature[local] = expected_temperature(gids[local]);

      for (std::size_t level = 0; level < nlevels; ++level) {
        moisture[local * nlevels + level] =
            expected_moisture(gids[local], level);
      }
    }
  }

  // Explicitly initialize destination tags too, so stale/default values
  // cannot accidentally satisfy the checks.
  adapter.write_tag("temperature", temperature);
  adapter.write_tag("soil_moisture", moisture);

  // SANITY CHECK
  std::vector<double> source_check(moisture.size(), sentinel);
  adapter.read_tag("soil_moisture", source_check);
  CHECK(source_check == moisture);

  if (is_forcing) {
    transport.send_exports("forcing");
    transport.complete_sends("forcing");
  } else {
    transport.receive_imports("land");
  }
  adapter.read_tag("temperature", temperature);
  adapter.read_tag("soil_moisture", moisture);

  // Finish communication on every rank before any value assertions.
  // This barrier is test scaffolding, not required by Transport.
  MPI_Barrier(MPI_COMM_WORLD);

  if (!is_forcing) {
    for (std::size_t local = 0; local < nlocal; ++local) {
      INFO("Global column ID: " << gids[local]);

      CHECK(temperature[local] == expected_temperature(gids[local]));

      for (std::size_t level = 0; level < nlevels; ++level) {
        INFO("Level: " << level);

        CHECK(moisture[local * nlevels + level] ==
              expected_moisture(gids[local], level));
      }
    }
  }
}
