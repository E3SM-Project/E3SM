/**
 * @file: hook.hpp
 *
 */
#ifndef E3SM_EMULATOR_HOOK_HPP
#define E3SM_EMULATOR_HOOK_HPP

#include "external_process.hpp"
#include <field_registry.hpp>

#include <emulator_config.hpp>
#include <field_types.hpp>
#include <grid_types.hpp>
#include <process.hpp>

#include <cstdint>
#include <filesystem>
#include <functional>
#include <span>
#include <string>
#include <variant>
#include <vector>

namespace e3sm::emulator {
using namespace coupler;

/**
  // ProcessRunOpts run_opts
  // GridDesc grid
  // EmulatorConfig emulator_config
*/
struct ProcessDescriptor {
  ProcessRunOpts run_opts;
  GridDesc grid;
  EmulatorConfig emulator_config;
};

// info needed to add a field to the hook
//   name - field name
//   Role - external binding field to look into.
struct HookedFieldDesc {
  std::string name;
  FieldRole role;
};

using MaskType = std::vector<std::uint8_t>;
using CallbackType = std::function<void(MaskType)>;

struct HookedField {

  HookedFieldDesc desc;
  FieldView data;

  // mask is vector instead of span to allow temporary masks
  // but then there needs to be a way to update them via the api.
  // We could force the external process to ensure the mask has
  // the same lifetime as the hook.
  std::vector<std::uint8_t> mask;

  // API decision:  for runtime dynamic masks (or other use cases)
  // we should make these functions take a HookContext struct as an argument.
  // but only if the needed is more common than not.
  std::function<void(MaskType)> pack;
  std::function<void(MaskType)> unpack;
};

/**
 * Hook:
 * - std::string name;
 * - ProcessDescriptor process_desc;
 * - std::vector<HookedField> fields;
 */
struct Hook {
  std::string name;
  ProcessDescriptor process_desc;
  std::vector<HookedField> fields;

  void add_field(HookedFieldDesc desc, FieldView data, MaskType mask,
                 CallbackType pack, CallbackType unpack);

  // Convenient overloads:
  void add_field(HookedFieldDesc desc, FieldView data) {
    add_field(
        std::move(desc), data, {}, [](MaskType mask) {}, [](MaskType mask) {});
  }
};

} // namespace e3sm::emulator
#endif
