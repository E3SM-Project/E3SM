#include "hook.hpp"

#include <field_registry.hpp>
#include <field_types.hpp>

#include <span>

namespace e3sm::emulator {

void Hook::add_field(HookedFieldDesc desc, FieldView data, MaskType mask,
                     CallbackType pack, CallbackType unpack) {

  fields.push_back({.desc = std::move(desc),
                    .data = data,
                    .mask = std::move(mask),
                    .pack = std::move(pack),
                    .unpack = std::move(unpack)});
}
} // namespace e3sm::emulator
