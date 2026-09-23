#include <field_registry.hpp>
#include <sstream>
#include <stdexcept>
#include <string_view>

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

  const auto pair = lookup_.find(RegistryKey{component, field_name});
  if (pair == lookup_.end()) {
    throw std::runtime_error("Coulding find " + field_name +
                             " in registry for " + component);
  }
  return pair->second;
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

const RegisteredFieldAttributes
read_attributes(const ekat::ParameterList& params) {
  return RegisteredFieldAttributes{
      .name = params.get<std::string>("attname"),
      .long_name = params.get<std::string>("longname"),
      .standard_name = params.get<std::string>("stdname"),
      .units = params.get<std::string>("units")};
}

std::string to_string(const RegisteredFieldAttributes& attr,
                      std::string_view spaces) {
  std::ostringstream out_str;
  out_str << spaces << "attname: " + attr.name + '\n'
          << spaces << "longname: " + attr.long_name + '\n'
          << spaces << "stdname: " + attr.standard_name + '\n'
          << spaces << "units: " + attr.units + '\n';
  return out_str.str();
}

} // namespace e3sm::coupler
