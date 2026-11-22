#ifndef E3SM_COUPLER_TESTS_FAKE_COMPONENTS_HPP
#define E3SM_COUPLER_TESTS_FAKE_COMPONENTS_HPP

#include "coupler.hpp"
#include "field_types.hpp"
#include <field_registry.hpp>
#include <process.hpp>
#include <string_view>
#include <vector>

namespace e3sm::emulator::coupler::test {

class FakeAtmosphere {
public:
  explicit FakeAtmosphere(std::size_t ncols)
      : temperature_(ncols, 0.0), precipitation_(ncols, 0.0),
        surface_flux_(ncols, 0.0) {}

  static constexpr std::string_view name() { return "atm"; }

  void populate_registry(FieldRegistry& registry) {

    coupling_.register_coupled_field(registry,
                                     RegisteredField{
                                         .role = FieldRole::Export,
                                         .component = "atm",
                                         .name = "temperature",
                                         .size = temperature_.size(),
                                     },
                                     FieldView{temperature_});

    coupling_.register_coupled_field(registry,
                                     {
                                         .role = FieldRole::Export,
                                         .component = "atm",
                                         .name = "precipitation",
                                         .size = precipitation_.size(),
                                     },
                                     FieldView{precipitation_});

    coupling_.register_coupled_field(registry,
                                     {
                                         .role = FieldRole::Import,
                                         .component = "atm",
                                         .name = "surface_flux",
                                         .size = surface_flux_.size(),
                                     },
                                     FieldView{surface_flux_});
  }

  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }

  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }

  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }

  void run(const ProcessRunOpts opts) {
    for (std::size_t i = 0; i < temperature_.size(); ++i) {
      temperature_[i] = 280.0 + static_cast<double>(i);
      precipitation_[i] = 0.1 * static_cast<double>(i);
    }
  }
  std::span<const int> get_local_gids() { return coupling_.get_local_gids(); }

  const std::vector<double>& temperature() const noexcept {
    return temperature_;
  }

  const std::vector<double>& precipitation() const noexcept {
    return precipitation_;
  }

  const std::vector<double>& surface_flux() const noexcept {
    return surface_flux_;
  }

private:
  std::vector<double> temperature_;
  std::vector<double> precipitation_;
  std::vector<double> surface_flux_;
  std::vector<FieldBinding> coupled_fields_;
  CouplingAdapter coupling_;
};

class FakeLand {
public:
  explicit FakeLand(std::size_t ncols)
      : temperature_(ncols, 0.0), precipitation_(ncols, 0.0),
        surface_flux_(ncols, 0.0) {}

  static constexpr std::string_view name() { return "lnd"; }
  void populate_registry(FieldRegistry& registry) {
    coupling_.register_coupled_field(registry,
                                     {
                                         .role = FieldRole::Import,
                                         .component = "lnd",
                                         .name = "temperature",
                                         .size = temperature_.size(),
                                     },
                                     FieldView{std::span{temperature_}});

    coupling_.register_coupled_field(registry,
                                     {
                                         .role = FieldRole::Import,
                                         .component = "lnd",
                                         .name = "precipitation",
                                         .size = precipitation_.size(),
                                     },
                                     FieldView{std::span{precipitation_}});

    coupling_.register_coupled_field(registry,
                                     {
                                         .role = FieldRole::Export,
                                         .component = "lnd",
                                         .name = "surface_flux",
                                         .size = surface_flux_.size(),
                                     },
                                     FieldView{std::span{surface_flux_}});
  }

  void run(const ProcessRunOpts) {
    for (std::size_t i = 0; i < surface_flux_.size(); ++i) {
      surface_flux_[i] = 2.0 * precipitation_[i];
    }
  }

  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }

  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }

  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() { return coupling_.get_local_gids(); }

  const std::vector<double>& temperature() const noexcept {
    return temperature_;
  }

  const std::vector<double>& precipitation() const noexcept {
    return precipitation_;
  }

  const std::vector<double>& surface_flux() const noexcept {
    return surface_flux_;
  }

private:
  std::vector<double> temperature_;
  std::vector<double> precipitation_;
  std::vector<double> surface_flux_;
  std::vector<FieldBinding> coupled_fields_;
  CouplingAdapter coupling_;
};

// One input and one output, each a single value. `model` computes the output
// from the input and how many times the component has run, from 1.
class RateComponent {
public:
  RateComponent(std::string name, std::string input, std::string output,
                std::function<double(double, int)> model)
      : name_(std::move(name)), input_name_(std::move(input)),
        output_name_(std::move(output)), model_(std::move(model)) {}

  std::string_view name() const { return name_; }
  void populate_registry(FieldRegistry& registry) {
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Import,
                                      .component = name_,
                                      .name = input_name_,
                                      .size = 1},
                                     input_);
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Export,
                                      .component = name_,
                                      .name = output_name_,
                                      .size = 1},
                                     output_);
  }
  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }
  void run(const ProcessRunOpts) { output_[0] = model_(input_[0], ++runs_); }
  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }
  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return coupling_.get_local_gids();
  }
  double input() const { return input_[0]; }

private:
  std::string name_;
  std::string input_name_;
  std::string output_name_;
  std::function<double(double, int)> model_;
  std::vector<double> input_{0.0};
  std::vector<double> output_{0.0};
  int runs_ = 0;
  CouplingAdapter coupling_;
};


// Binds one caller-owned span as `state` in both directions
class ExternalComponent {
public:
  explicit ExternalComponent(std::span<double> state) : state_(state) {}

  static constexpr std::string_view name() { return "external"; }
  void populate_registry(FieldRegistry& registry) {
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Import,
                                      .component = "external",
                                      .name = "state",
                                      .size = state_.size()},
                                     FieldView{state_});
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Export,
                                      .component = "external",
                                      .name = "state",
                                      .size = state_.size()},
                                     FieldView{state_});
  }
  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }
  void run(const ProcessRunOpts) {}
  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }
  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return coupling_.get_local_gids();
  }

private:
  std::span<double> state_;
  CouplingAdapter coupling_;
};

// Has `state` as both an input and an output, and advances it by one each run
class StatefulComponent {
public:
  explicit StatefulComponent(std::size_t ncols)
      : input_(ncols, 0.0), output_(ncols, 0.0) {}

  static constexpr std::string_view name() { return "stateful"; }
  void populate_registry(FieldRegistry& registry) {
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Import,
                                      .component = "stateful",
                                      .name = "state",
                                      .size = input_.size()},
                                     input_);
    coupling_.register_coupled_field(registry,
                                     {.role = FieldRole::Export,
                                      .component = "stateful",
                                      .name = "state",
                                      .size = output_.size()},
                                     output_);
  }
  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }
  void run(const ProcessRunOpts) {
    for (std::size_t i = 0; i < input_.size(); ++i) {
      output_[i] = input_[i] + 1.0;
    }
  }
  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }
  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return coupling_.get_local_gids();
  }

private:
  std::vector<double> input_;
  std::vector<double> output_;
  CouplingAdapter coupling_;
};


// Registers nothing, so it can be added to a driver without being in the plan
class UnplannedComponent {
public:
  static constexpr std::string_view name() { return "unplanned"; }
  void populate_registry(FieldRegistry&) {}
  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }
  void run(const ProcessRunOpts) {}
  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }
  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return coupling_.get_local_gids();
  }

private:
  CouplingAdapter coupling_;
};

} // namespace e3sm::emulator::coupler::test

#endif
