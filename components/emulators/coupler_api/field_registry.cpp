#include <field_registry.hpp>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace e3sm::emulator::coupler {

/**
 * @brief: Registers available fields with the registry.
 */
FieldID FieldRegistry::register_field(RegisteredField field) {

  if (field.component.empty()) {
    throw std::invalid_argument("Cannot register field '" + field.name +
                                "' without a component name");
  }

  if (field.name.empty()) {
    throw std::invalid_argument("Component '" + field.component +
                                "' tried to register a field without a name");
  }

  const auto id = static_cast<FieldID>(fields_.size());
  RegistryKey key{field.component, field.name, field.role};
  if (lookup_.contains(key)) {
    throw std::runtime_error("Component '" + field.component +
                             "' registered field '" + field.name + "' as " +
                             role_to_string(field.role) + " twice");
  }

  lookup_.emplace(std::move(key), id);
  fields_.emplace_back(std::move(field));

  return id;
}

const RegisteredField& FieldRegistry::get(const std::string& component,
                                          const std::string& field_name,
                                          FieldRole role) const {
  return fields_[get_id(component, field_name, role)];
}

const RegisteredField& FieldRegistry::get(FieldID id) const {
  if (id >= fields_.size()) {
    throw std::out_of_range("Unknown field ID " + std::to_string(id) + "; " +
                            std::to_string(fields_.size()) +
                            " field(s) are registered");
  }
  return fields_[id];
}

FieldID FieldRegistry::get_id(const std::string& component,
                              const std::string& field_name,
                              FieldRole role) const {

  const auto pair = lookup_.find(RegistryKey{component, field_name, role});
  if (pair == lookup_.end()) {
    std::string registered;
    for (const auto& field : fields_) {
      if (field.component == component) {
        registered += (registered.empty() ? "" : ", ") + field.name + " (" +
                      role_to_string(field.role) + ")";
      }
    }
    if (registered.empty()) {
      throw std::runtime_error("Field '" + field_name +
                               "' references component '" + component +
                               "', which has registered no fields. Was it "
                               "added to the driver?");
    }
    throw std::runtime_error("Component '" + component +
                             "' has not registered field '" + field_name +
                             "' as " + role_to_string(role) +
                             ". Registered fields: " + registered);
  }
  return pair->second;
}

bool FieldRegistry::contains(const std::string& component,
                             const std::string& field_name,
                             FieldRole role) const {
  RegistryKey key{component, field_name, role};
  return lookup_.contains(key);
}

std::string to_string(const MergeType merge_type) {
  switch (merge_type) {
  case MergeType::Direct:
    return "Direct";
  case MergeType::ScaledByFraction:
    return "ScaledByFraction";
  }
  throw std::runtime_error("Unsupported MergeType");
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

std::string role_to_string(const FieldRole role) {
  switch (role) {
  case FieldRole::Import:
    return "Import";
  case FieldRole::Export:
    return "Export";
  }
  throw std::runtime_error("Unsupported FieldRole");
}

std::string FieldRegistry::to_string() const {
  std::ostringstream out;
  out << "FieldRegistry: " << fields_.size() << " field(s)\n";

  for (FieldID id = 0; id < fields_.size(); ++id) {
    const auto& field = fields_[id];

    out << "  Field ID: " << id << '\n'
        << "    component: " << field.component << '\n'
        << "    role: " << role_to_string(field.role) << '\n'
        << "    size: " << field.size << '\n'
        << "    name: " << field.name << '\n';
  }

  return out.str();
}

} // namespace e3sm::emulator::coupler
