#ifndef E3SM_EMULATOR_HOOK_HPP
#define E3SM_EMULATOR_HOOK_HPP

#include "external_process.hpp"
#include "harness.hpp"

#include <cstdint>
#include <filesystem>
#include <functional>
#include <span>
#include <string>
#include <variant>
#include <vector>

namespace e3sm::emulator {

using RealField = std::span<double>;
using IntField = std::span<int>;
struct HookedField {

  std::string field_name;

  // mask is vector instead of span to allow temporary masks
  // but then there needs to be a way to update them via the api.
  // We could force the external process to ensure the mask has
  // the same lifetime as the hook.
  std::vector<std::uint8_t> mask;

  // using a variant here to make std::vector<HookedField> more ergonomic
  // std::variant<RealField, IntField> data;

  ExternalBindings bindings;

  // API decision:  for runtime dynamic masks (or other use cases)
  // we should make these functions take a HookContext struct as an argument.
  // but only if the needed is more common than not.
  std::function<void()> pack;
  std::function<void()> unpack;
};

struct Hook {
  std::string name;
  ProcessDescriptor process;
  std::filesystem::path emulator_config_file;
  std::filesystem::path coupling_config_file;
  std::vector<HookedField> fields;
};

} // namespace e3sm::emulator
#endif
