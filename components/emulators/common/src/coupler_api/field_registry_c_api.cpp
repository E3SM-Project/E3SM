
#include "coupler_types.hpp"
#include "field_registry.hpp"
#include <stdexcept>
#include <string>

namespace {
e3sm::coupler::FieldRole set_role(std::string role) {
  // will likey make this a fortran/c enum
  if (role == "source") {
    return e3sm::coupler::FieldRole::Export;
  } else if (role == "destination") {
    return e3sm::coupler::FieldRole::Import;
  } else {
    throw std::invalid_argument("Incorrect Role specified for Coupler Field");
  }
}

} // namespace

extern "C" {

void register_field(void* handle, const RegisteredFieldDesc* desc) {

  auto* registry = static_cast<e3sm::coupler::FieldRegistry*>(handle);
  e3sm::coupler::RegisteredField field{
      .role = set_role(desc->role),
      .attributes =
          {
              .component = desc->component,
              .name = desc->attributes.name,
              .long_name = desc->attributes.long_name,
              .standard_name = desc->attributes.standard_name,
              .units = desc->attributes.units,
          },
      .data = desc->data,
  };

  registry->register_field(std::move(field));
}


}

