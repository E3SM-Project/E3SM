#ifndef E3SM_COUPLER_API_FIELD_REGISTRY_HPP
#define E3SM_COUPLER_API_FIELD_REGISTRY_HPP
#include <string>
#include <unordered_map>
#include <vector>
namespace e3sm::coupler {

enum class FieldRole {
  Import,
  Export,
  // is this needed? there are fields that every component sends to the coupler
  // which merges them, and then each component has the merged state available
  // next timestep?
  ImportExport,
};

enum class MergeType {
  Direct,
  ScaledByFraction,
};

/**
 * @brief (C++ version) Attributes of Field to be registered
 * Fields:
 * - component: component that defines it
 * - name: name of field
 * - long_name: long name for field
 * - standard_name: standardized name for field
 * - units: units
 */
struct RegisteredFieldAttributes {
  std::string name;
  std::string long_name;
  std::string standard_name;
  std::string units;
};

struct RegisteredField {
  FieldRole role;
  std::string component;
  RegisteredFieldAttributes attributes;
  // TODO: Figure out how to register or couple fields that are downscaled or
  // otherwise transformed

  // TODO: Does this need to be generic or do we only have doubles?
  double* data = nullptr;
  std::size_t size = 0;
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
struct CoupledFieldEntry{
  std::string id;
  MergeType merge_type;
  RegisteredFieldAttributes attributes;
  std::vector<std::string> sources;
  std::vector<std::string> destinations;
};

/**
 * @brief Description of how a variable is coupled between components
 * Fields:
 *  - merge_type:  Enum for how multiple sources are merged
 *  - sources: fields that contribute to the same coupled state/flux
 *  - destinations: fields that consume the coupled state/flux
 */
struct CouplingRoute {
  MergeType merge_type;
  std::vector<const RegisteredField*> sources;
  std::vector<const RegisteredField*> destinations;
};

/**
 * @brief Registry that tells the coupler what fields are available to be
 * coupled Fields:
 * - fields_: unordered map with key = (component name, field name) to each
 * RegisteredField
 */
class FieldRegistry {
public:
  void register_field(RegisteredField field);
  const RegisteredField& get(const std::string& component,
                             const std::string& field_name) const;
  bool contains(const std::string& component,
                const std::string& field_name) const;

private:
  struct RegistryKey {
    std::string component;
    std::string field_name;

    bool operator==(const RegistryKey& other) const {
      return component == other.component && field_name == other.field_name;
    }
  };
  // taken from cpp reference
  struct RegistryKeyHash {
    std::size_t operator()(const RegistryKey& key) const noexcept {
      std::size_t h1 = std::hash<std::string>{}(key.component);
      std::size_t h2 = std::hash<std::string>{}(key.field_name);
      return h1 ^ (h2 << 1); // or use boost::hash_combine
    }
  };

  std::unordered_map<RegistryKey, RegisteredField, RegistryKeyHash> fields_;
};

} // namespace e3sm::coupler
#endif
