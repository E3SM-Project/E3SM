#ifndef E3SM_COUPLER_HPP
#define E3SM_COUPLER_HPP

#include <components.hpp>
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



/**
 * @brief Description of how a variable is coupled between components
 * Fields:
 *  - merge_type:  Enum for how multiple sources are merged
 *  - sources: fields that contribute to the same coupled state/flux
 *  - destinations: fields that consume the coupled state/flux
 */
struct CouplingRoute {
  MergeType merge_type;
  std::vector<FieldID> sources;
  std::vector<FieldID> destinations;
};


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

struct ActiveCouplingFields{
  std::vector<ImportBuffer> import_fields;
  std::vector<ExportBuffer> export_fields;
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
  const ActiveCouplingFields& coupling_plan(const std::string& component_name) const;

private:
  FieldRegistry registry_;
  std::vector<CouplingRoute> routes_;
  std::unordered_map<std::string, ActiveCouplingFields> coupling_plans_;

  std::vector<CoupledFieldEntry>
  read_coupling_fields_from_yaml(const std::filesystem::path& filename);
};


} // namespace e3sm::coupler
#endif
