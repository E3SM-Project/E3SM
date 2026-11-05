#ifndef E3SM_EMULATOR_MOAB_TRANSPORT_HPP
#define E3SM_EMULATOR_MOAB_TRANSPORT_HPP

#include "moab_adapter.hpp"
#include "mpi.h"
#include <span>
#include <string>
#include <vector>

namespace e3sm::coupler {

/*
  std::string field;
  std::string source_component;
  std::string destination_component;
*/
struct TransportRoute {
  std::string field;
  std::string source_component;
  std::string destination_component;
};

/*
 * int id -- app id
 * std::vector<int> ranks -- coupling_comm_ ranks
 */
struct ComponentInfo {
  int id = -1;
  std::vector<int> ranks;
};
  using ComponentMap = std::map<std::string, ComponentInfo>;
class Transport {
  using AdapterMap = std::map<std::string, MoabAdapter*>;

public:
  explicit Transport(MPI_Comm coupling_comm);
  ~Transport();
  Transport(const Transport&) = delete;
  Transport& operator=(const Transport&) = delete;
  Transport(Transport&&) = delete;
  Transport& operator=(Transport&&) = delete;

  void initialize(std::span<MoabAdapter* const> local_adapters,
                  std::span<const TransportRoute> routes);

  void send_exports(std::string_view component);
  void receive_imports(std::string_view component);
  void complete_sends(std::string_view component);
  ComponentMap discover_components(std::span<const TransportRoute> routes);

private:
  struct ConnectionKey {
    std::string source;
    std::string destination;

    bool operator==(const ConnectionKey& other) const {
      return source == other.source && destination == other.destination;
    }
  };
  // taken from cpp reference
  struct ConnectionKeyHash {
    std::size_t operator()(const ConnectionKey& key) const noexcept {
      std::size_t h1 = std::hash<std::string>{}(key.source);
      std::size_t h2 = std::hash<std::string>{}(key.destination);
      return h1 ^ (h2 << 1); // or use boost::hash_combine
    }
  };

  struct Connection {
    std::vector<std::string> fields;
    std::string tag_names; // colon separated list
    //  Mpi stuff
    int source_id = -1;
    int destination_id = -1;
    MPI_Group source_group = MPI_GROUP_NULL;
    MPI_Group destination_group = MPI_GROUP_NULL;

    bool send_pending = false;
  };

  void create_route_to_fields_map(std::span<const TransportRoute> routes);
  MPI_Group make_mpi_group(std::span<const int> ranks);

  MPI_Comm coupling_comm_ = MPI_COMM_NULL;
  std::map<std::string, MoabAdapter*> local_adapters_;

  using ConnectionMap =
      std::unordered_map<ConnectionKey, Connection, ConnectionKeyHash>;

  ConnectionMap connections_map_;
  bool initialized_ = false;
};

} // namespace e3sm::coupler
#endif
