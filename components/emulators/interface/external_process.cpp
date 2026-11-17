#include "external_process.hpp"

#include <field_registry.hpp>

#include <algorithm>
#include <iterator>
#include <stdexcept>

namespace e3sm::emulator {
void ExternalProcess::attach_fields(ExternalBindings bindings) {

  std::ranges::move(bindings.to_external,
                    std::back_inserter(bindings_.to_external));

  std::ranges::move(bindings.from_external,
                    std::back_inserter(bindings_.from_external));
}

void ExternalProcess::populate_registry(FieldRegistry& registry) {

  if (bindings_.from_external.empty() or bindings_.to_external.empty()) {
    throw std::runtime_error(
        "ExternalProcess: attach fields before populating the registry");
  }

  for (auto& field : bindings_.to_external) {
    coupling_.register_coupled_field(
        registry,
        RegisteredField{.role = FieldRole::Import,
                        .component = std::string(name()),
                        .name = field.name,
                        .size = field.data.size()},
        field.data);
  }

  for (auto& field : bindings_.from_external) {
    coupling_.register_coupled_field(
        registry,
        RegisteredField{.role = FieldRole::Export,
                        .component = std::string(name()),
                        .name = field.name,
                        .size = field.data.size()},
        field.data);
  }
}

void ExternalProcess::print_binding(std::ostream& os,
                                    const ExternalField& field) const {

  os << "      " << "::" << field.name << " ["
     << " size=" << field.data.size()
     << ", addr=" << static_cast<const void*>(field.data.data()) << "]\n";

  os << "       [ ";
  for (std::size_t i = 0; i < field.data.size(); i++) {
    os << field.data[i] << ", ";
  }
  os << "]\n";
}

void ExternalProcess::print_bindings(std::ostream& os) const {
  os << " ---- to_external ----- \n";
  for (const auto& binding : bindings_.to_external) {
    print_binding(os, binding);
  }

  os << " ----- from_external ------ \n";
  for (const auto& binding : bindings_.from_external) {
    print_binding(os, binding);
  }
  os << '\n';
}

} // namespace e3sm::emulator
