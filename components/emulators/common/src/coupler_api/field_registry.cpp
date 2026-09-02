#include <field_registry.hpp>
#include <stdexcept>

namespace e3sm::coupler {

/**
 * @brief: Registers available fields with the registry.
 */
FieldID FieldRegistry::register_field(RegisteredField field) {

  if (field.component.empty()) {
    throw std::invalid_argument(
        "Attempting to register field with no Component");
  }

  if (field.attributes.name.empty()) {
    throw std::invalid_argument("Attempting to register field with no name");
  }

  if (field.attributes.units.empty()) {
    throw std::invalid_argument("Attempting to register field with no units");
  }

  const auto id = static_cast<FieldID>(fields_.size());
  RegistryKey key{field.component, field.attributes.name};
  if (lookup_.contains(key)) {
    throw std::runtime_error("Attempted to register field twice");
  }

  auto [it, inserted] = lookup_.emplace(std::move(key), id);
  fields_.emplace_back(std::move(field));

  return id;
}

const RegisteredField& FieldRegistry::get(const std::string& component,
                                          const std::string& field_name) const {
  RegistryKey key{component, field_name};
  auto id = lookup_.find(key)->second;
  return fields_.at(id);
}

const RegisteredField& FieldRegistry::get(FieldID id) const {
  return fields_.at(id);
}

FieldID FieldRegistry::get_id(const std::string& component,
                              const std::string& field_name) const {
  return lookup_.find(RegistryKey{component, field_name})->second;
}

bool FieldRegistry::contains(const std::string& component,
                             const std::string& field_name) const {
  RegistryKey key{component, field_name};
  return lookup_.contains(key);
}

std::string to_string(const MergeType merge_type) {
  switch (merge_type) {
  case MergeType::Direct:
    return "Direct";
  case MergeType::ScaledByFraction:
    return "ScaledByFraction";
  }
}
} // namespace e3sm::coupler
