#ifndef E3SM_COUPLER_HPP
#define E3SM_COUPLER_HPP

#include "field_registry.hpp"
#include "field_types.hpp"

#include <filesystem>
#include <ostream>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>

namespace e3sm::emulator::coupler {

/**
 * @brief Holds yaml entries for active coupled fields
 * Fields:
 * - string id: unique label for the route
 * - string field_name: name of the field the route carries
 * - MergeType merge_type
 * - vector<string> sources: list of component sources
 * - vector<string> destinations: list of component destinations
 */
struct CoupledFieldEntry {
  std::string id;
  std::string field_name;
  MergeType merge_type;
  std::vector<std::string> sources;
  std::vector<std::string> destinations;
};

std::vector<CoupledFieldEntry>
read_coupling_fields_from_yaml(const std::filesystem::path& filename);

std::string to_string(const CoupledFieldEntry& entry);

inline std::ostream& operator<<(std::ostream& os,
                                const CoupledFieldEntry& entry) {
  return os << to_string(entry) << '\n';
}

/**
 * @brief Description of how a variable is coupled between components
 * Fields:
 *  - merge_type:  Enum for how multiple sources are merged
 *  - std::string name: unique label for the route
 *  - std::string field_name: name of the field the route carries
 *  - sources: fields that contribute to the same coupled state/flux
 *  - destinations: fields that consume the coupled state/flux
 *  - owned_buffer: Optional Owned data buffer that components copy to/from
 *  - buffer: span pointing to the buffer to be used by the route
 *  - source_component, destination_components: names of the components involved
 *  - averages: destinations run less often than the source, so they import the
 *    mean of the source's exports since they last ran, rather than the latest
 *  - accumulator, n_samples, mean: storage for that mean
 */
struct CouplingRoute {
  MergeType merge_type;
  std::string name;
  std::string field_name;
  std::vector<FieldID> sources; 
  std::vector<FieldID> destinations;


  FieldStorage owned_buffer;
  FieldView buffer;

  std::string source_component;
  std::vector<std::string> destination_components;

  // Accmulator details
  bool averages = false;
  std::size_t n_samples = 0;
  FieldStorage accumulator;
  FieldStorage mean;
};

/**
 * @brief Provides list of fields that are active
 * import_ids: list of FieldIDs that correspond to imports
 * export_ids: list of FieldIDs that correspond to exports
 */
struct ActiveCouplingFields {
  std::vector<FieldID> import_ids;
  std::vector<FieldID> export_ids;
};

/**
 * @brief Non-owning mapping to export and import buffers per component
 * fields:
  std::vector<ExportBuffer> exports;
  std::vector<ImportBuffer> imports;
 */
struct ComponentBuffers {
  std::vector<ExportBuffer> exports;
  std::vector<ImportBuffer> imports;
};

/**
 * @brief Coupler handles the exchange of variables between CoupledProcesses
 * Fields:
 * - registry_: main registry of fields that each component exposes for import
 * or export
 * - routes_: The active routes of source fields to destination fields
 * Functions:
 * - build_routes: create routes from CoupledFieldEntry list
 * config file
 */
using RouteBufferMap = std::unordered_map<std::string, FieldView>;

/// How many base steps pass between runs of each component. A component that
/// isn't listed, such as an external one, doesn't take part in averaging.
using ComponentIntervals = std::unordered_map<std::string, int>;

class Coupler {
public:
  Coupler() = default;

  // neither copyable nor movable
  Coupler(Coupler&) = delete;
  Coupler(const Coupler&) = delete;
  Coupler& operator=(const Coupler&) = delete;
  Coupler& operator=(const Coupler&&) = delete;

  FieldRegistry& registry() noexcept { return registry_; }
  const FieldRegistry& registry() const noexcept { return registry_; }

  // methods to construct routes and attach buffers used for import/export
  void build_routes(const std::vector<CoupledFieldEntry> &coupling_fields);
  void build_buffers(RouteBufferMap external_buffers,
                     const ComponentIntervals& intervals = {});

  /// Adds the component's latest exports to the means of its averaging routes.
  void accumulate_exports(std::string_view component_name);

  /// Computes the means that the component is about to import, and resets them.
  /// @throws std::runtime_error if nothing was exported since the last import
  void prepare_imports(std::string_view component_name);

  const ActiveCouplingFields&
  coupling_plan(const std::string& component_name) const;

  const std::span<ExportBuffer>
  export_buffers(std::string_view component_name) {
    return component_buffers(component_name).exports;
  }

  const std::span<ImportBuffer>
  import_buffers(std::string_view component_name) {
    return component_buffers(component_name).imports;
  }

  void print_coupling_plan(std::ostream& os) const;

  void print_field_buffer(std::ostream& os, FieldID id,
                          const ConstFieldView buffer,
                          std::size_t count = 0) const;

private:
  FieldRegistry registry_;
  std::vector<CouplingRoute> routes_;
  std::unordered_map<std::string, ActiveCouplingFields> coupling_plans_;
  std::unordered_map<std::string, ComponentBuffers> buffers_;

  ComponentBuffers& component_buffers(std::string_view component_name);
};



} // namespace e3sm::emulator::coupler
#endif
