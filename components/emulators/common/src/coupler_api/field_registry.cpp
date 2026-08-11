#include <field_registry.hpp>
#include <stdexcept>

namespace e3sm::coupler {

void FieldRegistry::register_field(RegisteredField field) {

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

  if (field.data == nullptr) {
    throw std::invalid_argument("Attempting to Register field with nullptr");
  }

  RegistryKey key{field.component, field.attributes.name};
  auto [it, inserted] = fields_.emplace(std::move(key), std::move(field));

  if (!inserted) {
    // throw for now unless we want to handle potentially stale pointers?
    throw std::runtime_error("Attempted to register field twice");
  }
  return;
}

const RegisteredField& FieldRegistry::get(const std::string &component,
                                          const std::string &field_name) const {
  RegistryKey key{component, field_name};
  return fields_.find(key)->second;
}

bool FieldRegistry::contains(const std::string &component,
                             const std::string &field_name) const {
  RegistryKey key{component, field_name};
  return fields_.contains(key);
}

} // namespace e3sm::coupler
