#ifndef E3SM_EMULATOR_CPL7_HARNESS_HPP
#define E3SM_EMULATOR_CPL7_HARNESS_HPP

#include "coupler_driver.hpp"
#include "emulator.hpp"
#include "emulator_config.hpp"
#include "grid_types.hpp"
#include "process.hpp"
#include <coupler.hpp>
#include <cpl7_component.hpp>
#include <map>

namespace e3sm::emulator {

using namespace e3sm::coupler;

/**
  // ProcessRunOpts run_opts;
  // GridDesc grid;
  // EmulatorConfig emulator_config;
  // 
*/
struct ProcessDescriptor {
  ProcessRunOpts run_opts;
  GridDesc grid;
  EmulatorConfig emulator_config;
};

class Harness {
  using ProcessMap = std::map<std::string, ProcessRunOpts>;

public:
  // neither copyable nor movable
  Harness() = default;
  Harness(Harness&) = delete;
  Harness(const Harness&) = delete;
  Harness& operator=(const Harness&) = delete;
  Harness& operator=(const Harness&&) = delete;

  void transfer_bindings(ExternalBindings bindings);

  void initialize(const std::filesystem::path& plan);
  void run_process(std::string_view name);
  void add_process(ProcessDescriptor desc);

private:
  ProcessMap processes_;
  CouplerDriver driver_;
  CPL7Component* external_component_ = nullptr;
  bool initialized_ = false;
};

} // namespace e3sm::emulator

#endif
