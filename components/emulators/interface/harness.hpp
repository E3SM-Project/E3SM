#ifndef E3SM_EMULATOR_CPL7_HARNESS_HPP
#define E3SM_EMULATOR_CPL7_HARNESS_HPP

/**
 * @file: harness.hpp
 * @brief: Implements `Harness` class for bridging external processes to processes
 * owned by the CouplerDriver.
 *
 * Design:
 *  ____________          ___________         ___________
 * |            |        |           |       |          |
 * |  External  |<-----> |  Harness  |<----->|  Coupler |
 * |  Processes |        |___________|       |  Driver  |
 * |____________|                            |__________|
 *
 * The Harness contains a pointer to an external process and owns a
 * CouplerDriver instance.
 *
 * `to_external` --> run the import fields method on the external process
 * `from_external` --> run the export fields method on the external process
 *
 * Memory buffers from the external process are passed through to CouplerDriver
 * via `ExternalBindings`
 *
 * Processes are added via a `ProcessDescriptor`, which holds the information
 * for configuring and running the process
 */

#include <coupler.hpp>
#include <coupler_driver.hpp>
#include <cpl7_component.hpp>
#include <emulator_config.hpp>
#include <grid_types.hpp>
#include <process.hpp>

namespace e3sm::emulator {

using namespace coupler;

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
  struct ProcessState {
    ProcessRunOpts run_opts;
    int calls_since_run = 0;
  };
  using ProcessMap = std::map<std::string, ProcessState>;

public:
  // neither copyable nor movable
  Harness() = default;
  Harness(Harness&) = delete;
  Harness(const Harness&) = delete;
  Harness& operator=(const Harness&) = delete;
  Harness& operator=(const Harness&&) = delete;

  void transfer_bindings(ExternalBindings bindings);

  void initialize(const std::filesystem::path& plan, RouteBufferMap external);
  void initialize(const std::filesystem::path& plan) {
    initialize(plan, RouteBufferMap{});
  }
  /// Call once per base step. A process with a coupling_interval of N imports,
  /// runs, and exports on every Nth call.
  /// @return whether the process ran
  bool run_process(std::string_view name);
  void add_process(ProcessDescriptor desc);
  void import_to_external(std::string_view name);
  void export_from_external(std::string_view name);

  /// Diagnostic: writes every value of every external binding to os.
  void print_bindings(std::ostream& os) const;

  template <CoupledProcess C> C& get_process(std::string_view name) {
    return driver_.get_component_as<C>(name);
  };

private:
  ProcessMap processes_;
  CouplerDriver driver_;
  CPL7Component* external_component_ = nullptr;
  bool initialized_ = false;
};

} // namespace e3sm::emulator

#endif
