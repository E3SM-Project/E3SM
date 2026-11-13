#include "moab_transport.hpp"

#include <moab/Types.hpp>
#include <moab/iMOAB.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <limits>
#include <ranges>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace e3sm::emulator::coupler {
namespace {

int checked_size_to_int(std::size_t size) {
  if (size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error("Size excedes int length");
  }
  return static_cast<int>(size);
}
// Fatal errors abort the joint communicator so peers aren't left waiting.
// Recoverable/configuration errors can gain a collective error path later.
void check_mpi(int err, MPI_Comm comm, const char* operation) {
  if (err != MPI_SUCCESS) {
    std::cerr << "MPI failure: " << operation << '\n';
    MPI_Abort(comm, err);
    std::terminate();
  }
}

void check_moab(ErrCode err, MPI_Comm comm, const char* operation) {
  if (err != moab::MB_SUCCESS) {
    std::cerr << "MOAB failure: " << operation
              << " code=" << static_cast<int>(err) << '\n';
    MPI_Abort(comm, static_cast<int>(err));
    std::terminate();
  }
}

std::string colon_delimited_string(std::span<const std::string> fields) {
  std::string result = "";
  if (fields.empty()) {
    return result;
  }
  result = fields.front();
  for (const auto& field : fields | std::views::drop(1)) {
    result += ":" + field;
  }
  return result;
}

} // namespace

Transport::Transport(MPI_Comm coupling_comm) : coupling_comm_(coupling_comm) {
  if (coupling_comm_ == MPI_COMM_NULL) {
    throw std::invalid_argument("Transport needs a valid communicator");
  }
}

void Transport::create_route_to_fields_map(
    std::span<const TransportRoute> routes) {

  for (const TransportRoute& route : routes) {
    const ConnectionKey key{.source = route.source_component,
                            .destination = route.destination_component};

    ConnectionMap::iterator it = connections_map_.find(key);
    if (it == connections_map_.end()) {
      Connection new_connection;
      new_connection.fields.push_back(route.field);
      connections_map_[key] = std::move(new_connection);
    } else {
      it->second.fields.push_back(route.field);
    }
  }
}

/*
 * Each component needs to know the moab app ids to exchange with
 * for moab to compute the communication graph
 * Let's assume we have three components:
 * - forcing and land on the same 2 ranks
 * - atm on 4 ranks
 *
 *   land+forcing:  0 1
 *   atm:  2 3 4 5
 *
 *   This function generates a num ranks X num components table to look up
 *   which ranks hold which components
 */
ComponentMap
Transport::discover_components(std::span<const TransportRoute> routes) {

  ComponentMap components;
  for (const TransportRoute& route : routes) {
    components.try_emplace(route.source_component);
    components.try_emplace(route.destination_component);
  }

  std::vector<int> local_component_ids;
  local_component_ids.reserve(components.size());

  bool valid_local_ids = true;

  for (const auto& [name, info] : components) {
    const auto it = local_adapters_.find(name);
    const int id = it == local_adapters_.end() ? 0 : it->second->component_id();

    if (it != local_adapters_.end() && id <= 0) {
      valid_local_ids = false;
    }

    local_component_ids.push_back(id);
  }
  int local_ok = valid_local_ids ? 1 : 0;
  int all_ok = 0;

  check_mpi(
      MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, coupling_comm_),
      coupling_comm_, "validating component IDs");

  if (!all_ok) {
    throw std::runtime_error("Component IDs must be positive");
  }

  int comm_size = 0;
  check_mpi(MPI_Comm_size(coupling_comm_, &comm_size), coupling_comm_,
            "getting communicator size");

  const int count = checked_size_to_int(local_component_ids.size());
  const std::size_t width = local_component_ids.size();

  // buffer to recieve complete list of ids in the comm
  //  This is our table.
  //  Using the example above we will populate it for ids: [atm, lnd, for] =
  //  [1,2,3]
  //
  //  Rank   | for | lnd | atm
  //  0      | 3   |  2  | 0
  //  1      | 3   | 2   | 0
  //  2,3,4,5|  0  |  0  | 1
  std::vector<int> gathered_component_ids(static_cast<std::size_t>(comm_size) *
                                          width);

  check_mpi(MPI_Allgather(local_component_ids.data(), count, MPI_INT,
                          gathered_component_ids.data(), count, MPI_INT,
                          coupling_comm_),
            coupling_comm_, "gathering component membership");

  // Loop through each rank.
  for (std::size_t rank = 0; rank < comm_size; ++rank) {
    const auto rank_ids = std::span<const int>{gathered_component_ids}.subspan(
        rank * width, width);

    std::size_t column = 0;
    for (auto& [name, info] : components) {
      const int id = rank_ids[column++];

      if (id == 0) {
        continue;
      }

      if (info.id != -1 && info.id != id) {
        throw std::runtime_error("Inconsistent component ID for " + name);
      }

      info.id = id;
      info.ranks.push_back(rank);
    }
  }

  for (const auto& [name, info] : components) {
    if (info.ranks.empty()) {
      throw std::runtime_error("No ranks host component " + name);
    }
  }

  return components;
}

MPI_Group Transport::make_mpi_group(std::span<const int> ranks) {

  MPI_Group couple_group = MPI_GROUP_NULL;
  MPI_Group rank_group = MPI_GROUP_NULL;

  check_mpi(MPI_Comm_group(coupling_comm_, &couple_group), coupling_comm_,
            "getting couple group");

  check_mpi(MPI_Group_incl(couple_group, static_cast<int>(ranks.size()),
                           ranks.data(), &rank_group),
            coupling_comm_, "creating rank group");

  check_mpi(MPI_Group_free(&couple_group), coupling_comm_,
            "freeing temp group");

  return rank_group;
}

void Transport::initialize(std::span<MoabAdapter* const> local_adapters,
                           std::span<const TransportRoute> routes) {

  if (initialized_) {
    throw std::logic_error("Transport already initialized");
  }
  if (routes.empty()) {
    throw std::invalid_argument("Transport requires at least one route");
  }

  for (MoabAdapter* adapter : local_adapters) {
    if (adapter == nullptr) {
      throw std::invalid_argument("Null local adapter");
    }

    const auto [it, inserted] =
        local_adapters_.emplace(std::string(adapter->name()), adapter);

    if (!inserted) {
      throw std::invalid_argument("Duplicate local component");
    }
  }

  create_route_to_fields_map(routes);
  const ComponentMap components = discover_components(routes);

  // ComputeCommGraph's selector: 2 means match vertex GLOBAL_IDs.
  int source_type = 2;
  int destination_type = 2;

  using EndPointPair = std::pair<std::string, std::string>;
  std::set<EndPointPair> built_pairs;

  for (auto& [endpoints, connection] : connections_map_) {
    connection.tag_names = colon_delimited_string(connection.fields);
    const ComponentInfo& source = components.at(endpoints.source);
    const ComponentInfo& dest = components.at(endpoints.destination);

    connection.source_id = source.id;
    connection.destination_id = dest.id;
    connection.source_group = make_mpi_group(source.ranks);
    connection.destination_group = make_mpi_group(dest.ranks);

    const AdapterMap::const_iterator source_it =
        local_adapters_.find(endpoints.source);

    const AdapterMap::const_iterator dest_it =
        local_adapters_.find(endpoints.destination);

    int source_app =
        source_it == local_adapters_.end() ? -1 : source_it->second->app_id();

    int dest_app =
        dest_it == local_adapters_.end() ? -1 : dest_it->second->app_id();

    const EndPointPair pair{std::min(endpoints.source, endpoints.destination),
                            std::max(endpoints.source, endpoints.destination)};

    if (!built_pairs.contains(pair)) {

      check_moab(iMOAB_ComputeCommGraph(
                     &source_app, &dest_app, &coupling_comm_,
                     &connection.source_group, &connection.destination_group,
                     &source_type, &destination_type, &connection.source_id,
                     &connection.destination_id),
                 coupling_comm_, "computing communication graph");

      built_pairs.insert(pair);
    }
  }
  initialized_ = true;
}

void Transport::send_exports(std::string_view component) {

  if (!initialized_) {
    throw std::logic_error("Transport not initialized");
  }

  const AdapterMap::const_iterator adapter =
      local_adapters_.find(std::string(component));

  if (adapter == local_adapters_.end()) {
    return; // Component is not hosted on this rank.
  }

  for (auto& [endpoints, connection] : connections_map_) {
    if (endpoints.source != component) {
      continue;
    }

    if (connection.send_pending) {
      throw std::logic_error("Previous send has not been completed");
    }

    int app_id = adapter->second->app_id();
    int context = connection.destination_id;

    check_moab(iMOAB_SendElementTag(&app_id, connection.tag_names.c_str(),
                                    &coupling_comm_, &context),
               coupling_comm_, "sending field tags");

    connection.send_pending = true;
  }
}

void Transport::receive_imports(std::string_view component) {
  if (!initialized_) {
    throw std::logic_error("Transport not initialized");
  }

  const AdapterMap::const_iterator adapter =
      local_adapters_.find(std::string(component));

  if (adapter == local_adapters_.end()) {
    return;
  }

  for (auto& [endpoints, connection] : connections_map_) {
    if (endpoints.destination != component) {
      continue;
    }

    int app_id = adapter->second->app_id();
    int context = connection.source_id;

    check_moab(iMOAB_ReceiveElementTag(&app_id, connection.tag_names.c_str(),
                                       &coupling_comm_, &context),
               coupling_comm_, "receiving field tags");
  }
}

void Transport::complete_sends(std::string_view component) {
  if (!initialized_) {
    throw std::logic_error("Transport not initialized");
  }

  const AdapterMap::const_iterator adapter =
      local_adapters_.find(std::string(component));

  if (adapter == local_adapters_.end()) {
    return;
  }

  for (auto& [endpoints, connection] : connections_map_) {
    if (endpoints.source != component || !connection.send_pending) {
      continue;
    }

    int app_id = adapter->second->app_id();
    int context = connection.destination_id;

    check_moab(iMOAB_FreeSenderBuffers(&app_id, &context), coupling_comm_,
               "completing field sends");

    connection.send_pending = false;
  }
}

Transport::~Transport() {
  // Outstanding sends must already be completed.
  // Transport must be destroyed before MPI_Finalize().
  for (auto& [endpoints, connection] : connections_map_) {
    if (connection.source_group != MPI_GROUP_NULL) {
      MPI_Group_free(&connection.source_group);
    }

    if (connection.destination_group != MPI_GROUP_NULL) {
      MPI_Group_free(&connection.destination_group);
    }
  }
}

} // namespace e3sm::emulator::coupler
