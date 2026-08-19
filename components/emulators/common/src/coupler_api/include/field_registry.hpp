#ifndef E3SM_COUPLER_API_FIELD_REGISTRY_HPP
#define E3SM_COUPLER_API_FIELD_REGISTRY_HPP
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
namespace e3sm::coupler {

enum class FieldRole {
  Import,
  Export,
  ImportExport,
};

enum class MergeType {
  Direct,
  ScaledByFraction,
};

std::string to_string(const MergeType merge_type);

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

using FieldID = std::size_t;

template<typename T>
struct FieldBuffer{
  FieldID id;
  std::span<T> data;
};


using ImportBuffer = FieldBuffer<const double>;
using ExportBuffer = FieldBuffer<double>;

struct RegisteredField {
  FieldRole role;
  std::string component;
  RegisteredFieldAttributes attributes;
  std::size_t size = 0;
};

/**
 * @brief Registry that tells the coupler what fields are available to be
 * coupled Fields:
 * - fields_: unordered map with key = (component name, field name) to each
 * RegisteredField
 */
class FieldRegistry {
public:
  FieldRegistry() = default;

  FieldID register_field(RegisteredField field);
  const RegisteredField& get(const std::string& component,
                             const std::string& field_name) const;
  const RegisteredField& get(FieldID id) const;

  FieldID get_id(const std::string& component,
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

  std::vector<RegisteredField> fields_;
  std::unordered_map<RegistryKey, FieldID, RegistryKeyHash> lookup_;
};

} // namespace e3sm::coupler
#endif
