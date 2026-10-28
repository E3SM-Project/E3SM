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

struct BoundBuffers {
  std::span<double> out_buffer;
  std::span<double> in_buffer;

  std::vector<double> packed_out;
  std::vector<double> packed_in;

  std::vector<ExternalField> external_fields_;
};

struct ProcessDescriptor {
  ProcessRunOpts run_opts;
  GridDesc grid;
  EmulatorConfig emulator_config;
  std::vector<ExternalField> external_fields;
};

/**
 * ProcessDescriptor descriptors;
 * BoundBuffers buffers;
 */
struct BoundProcess {
  ProcessDescriptor descriptor;
  BoundBuffers buffer;
};

class Harness {
  using ProcessMap = std::map<std::string, BoundProcess>;

public:
  // neither copyable nor movable
  Harness() = default;
  Harness(Harness&) = delete;
  Harness(const Harness&) = delete;
  Harness& operator=(const Harness&) = delete;
  Harness& operator=(const Harness&&) = delete;

  void configure_emulator(EmulatorConfig config);
  void map_fields_to_buffer(std::vector<std::string> field_list);
  void run_process(std::string_view name);
  void add_process(ProcessDescriptor desc);

private:
  void component_send(std::string_view name);
  void component_receive(std::string_view name);
  void bind_buffers();

  ProcessMap processes_;
  CouplerDriver driver_;
  CPL7Component* external_component_ = nullptr;
};

} // namespace e3sm::emulator

#endif
