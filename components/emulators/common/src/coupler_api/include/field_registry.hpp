#ifndef COUPLER_API_FIELD_REGISTRY_HPP
#define COUPLER_API_FIELD_REGISTRY_HPP
#include <string>
#include <unordered_map>
#include <vector>
namespace e3sm::coupler {

enum class FieldRole {
  Source,
  Consumer,
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
  std::string component;
  std::string name;
  std::string long_name;
  std::string standard_name;
  std::string units;
};

struct RegisteredField {
  FieldRole role;
  RegisteredFieldAttributes attributes;

  // TODO: Does this need to be generic or do we only have doubles?
  double* data = nullptr;
  std::size_t size = 0;
};

struct CouplingRoute {
  MergeType merge_type;
  std::vector<const RegisteredField*> sources;
  std::vector<const RegisteredField*> destinations;
};

class FieldRegistry {
public:
  void register_field(RegisteredField field);


private:
  struct RegistryKey {
    std::string component;
    std::string field_name;

    // maybe this can be = default?
    bool operator==(const RegistryKey& other) const = default;
  };
  // taken from cpp reference
  struct MyHash {
    std::size_t operator()(const RegistryKey& key) const noexcept {
      std::size_t h1 = std::hash<std::string>{}(key.component);
      std::size_t h2 = std::hash<std::string>{}(key.field_name);
      return h1 ^ (h2 << 1); // or use boost::hash_combine
    }
  };

  std::unordered_map<RegistryKey, RegisteredField, MyHash> fields_;
};

} // namespace e3sm::coupler
#endif
