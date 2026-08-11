#ifndef E3SM_COUPLER_HPP
#define E3SM_COUPLER_HPP

#include <bits/std_thread.h>
#include <coupler_types.hpp>
#include <field_registry.hpp>
#include <ostream>
#include <sstream>
#include <string>

namespace e3sm::coupler {

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

class Coupler {
public:
  Coupler() = default;
  FieldRegistry& registry() noexcept { return registry_; }
  const FieldRegistry& registry() const noexcept { return registry_; }

  void build_routes();
  // void exchange();

  const std::vector<CoupledFieldEntry>
  read_coupling_fields_from_yaml(const std::string& filename);

private:
  FieldRegistry registry_;
  std::vector<CouplingRoute> routes_;
};

} // namespace e3sm::coupler
#endif
