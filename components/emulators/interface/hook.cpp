#include "hook.hpp"

#include <field_registry.hpp>

#include <span>

namespace e3sm::emulator {

void Hook::add_field(HookedFieldDesc desc, std::span<double> data,
                     MaskType mask, CallbackType pack, CallbackType unpack) {

  fields.push_back({.desc = std::move(desc),
                    .data = data,
                    .mask = std::move(mask),
                    .pack = std::move(pack),
                    .unpack = std::move(unpack)});
}
} // namespace e3sm::emulator
