#ifndef E3SM_COUPLER_HPP
#define E3SM_COUPLER_HPP

#include <concepts>
#include <coupler_types.hpp>
#include <field_registry.hpp>
#include <filesystem>
#include <ostream>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>

namespace e3sm::coupler {

using CouplingBuffer = FieldBuffer<double>;

/**
 * @brief Holds yaml entries for active coupled fields
 * Fields:
 * - merge_type
 * - id: unique label for coupled field
 * - attributes: Field metadata
 * - sources: list of component sources
 * - destinations: list of component destinations
 */
struct CoupledFieldEntry {
  std::string id;
  MergeType merge_type;
  RegisteredFieldAttributes attributes;
  std::vector<std::string> sources;
  std::vector<std::string> destinations;
};

std::string to_string(const CoupledFieldEntry& entry);

inline std::ostream& operator<<(std::ostream& os,
                                const CoupledFieldEntry& entry) {
  return os << to_string(entry) << '\n';
}

/**
 * @brief Description of how a variable is coupled between components
 * Fields:
 *  - merge_type:  Enum for how multiple sources are merged
 *  - sources: fields that contribute to the same coupled state/flux
 *  - destinations: fields that consume the coupled state/flux
 *  - buffers: Owned data buffers that components copy to/from
 */
struct CouplingRoute {
  MergeType merge_type;
  std::vector<FieldID> sources;
  std::vector<FieldID> destinations;
  std::vector<double> buffer;
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
 */
struct ComponentBuffers {
  std::vector<ExportBuffer> exports;
  std::vector<ImportBuffer> imports;
};

/**
 * @brief Coupler handles the exchange of variables between E3SMComponents
 * Fields:
 * - registry_: main registry of fields that each component exposes for import
 * or export
 * - routes_: The active routes of source fields to destination fields
 * Functions:
 * - build_routes: create routes from CoupledFieldEntry list
 * - read_coupling_fields_from_yaml: generate CoupledFieldEntry list from yaml
 * config file
 */
class Coupler {
public:
  Coupler() = default;
  FieldRegistry& registry() noexcept { return registry_; }
  const FieldRegistry& registry() const noexcept { return registry_; }
  void build_routes(const std::filesystem::path& filename);

  const ActiveCouplingFields&
  coupling_plan(const std::string& component_name) const;

  const std::span<ExportBuffer> export_buffers(std::string_view component_name) {
    return buffers_.at(std::string(component_name)).exports;
  }

  const std::span<ImportBuffer> import_buffers(std::string_view component_name) {
    return buffers_.at(std::string(component_name)).imports;
  }

  void print_coupling_plan(std::ostream& os) const;
  void print_field_buffer(std::ostream& os, FieldID id,
                          const std::span<const double> buffer, std::size_t count=0) const;

private:
  FieldRegistry registry_;
  std::vector<CouplingRoute> routes_;
  std::unordered_map<std::string, ActiveCouplingFields> coupling_plans_;
  std::unordered_map<std::string, ComponentBuffers> buffers_;

  std::vector<CoupledFieldEntry>
  read_coupling_fields_from_yaml(const std::filesystem::path& filename);

  void build_buffers();

};

} // namespace e3sm::coupler
#endif
