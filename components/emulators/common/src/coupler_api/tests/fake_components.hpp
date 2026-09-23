#ifndef E3SM_COUPLER_TESTS_FAKE_COMPONENTS_HPP
#define E3SM_COUPLER_TESTS_FAKE_COMPONENTS_HPP

#include "coupler.hpp"
#include <algorithm>
#include <field_registry.hpp>
#include <string_view>
#include <vector>
#include <process.hpp>

namespace e3sm::coupler::test {


class FakeAtmosphere {
public:
  explicit FakeAtmosphere(std::size_t ncols)
      : temperature_(ncols, 0.0), precipitation_(ncols, 0.0),
        surface_flux_(ncols, 0.0) {}

  static constexpr std::string_view name(){
    return "atm";
  }

  void populate_registry(FieldRegistry& registry) {

    coupling_.register_coupled_field(
        registry,
        RegisteredField{

            .role = FieldRole::Export,
            .component = "atm",
            .attributes =
                {
                    .name = "temperature",
                    .long_name = "Atmospheric temperature",
                    .standard_name = "air_temperature",
                    .units = "K",
                },
            .size = temperature_.size(),
        },
        temperature_);

    coupling_.register_coupled_field(
        registry,
        {
            .role = FieldRole::Export,
            .component = "atm",
            .attributes =
                {
                    .name = "precipitation",
                    .long_name = "Precipitation",
                    .standard_name = "precipitation_flux",
                    .units = "kg m-2 s-1",
                },
            .size = precipitation_.size(),
        },
        precipitation_);

    coupling_.register_coupled_field(
        registry,
        {
            .role = FieldRole::Import,
            .component = "atm",
            .attributes =
                {
                    .name = "surface_flux",
                    .long_name = "Surface energy flux",
                    .standard_name = "",
                    .units = "W m-2",
                },
            .size = surface_flux_.size(),
        },
        surface_flux_);
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

  static constexpr std::string_view name(){
    return "lnd";
  }
  void populate_registry(FieldRegistry& registry) {
    coupling_.register_coupled_field(
        registry,
        {
            .role = FieldRole::Import,
            .component = "lnd",
            .attributes =
                {
                    .name = "temperature",
                    .long_name = "Atmospheric temperature over land",
                    .standard_name = "air_temperature",
                    .units = "K",
                },
            .size = temperature_.size(),
        },
        temperature_);

    coupling_.register_coupled_field(
        registry,
        {
            .role = FieldRole::Import,
            .component = "lnd",
            .attributes =
                {
                    .name = "precipitation",
                    .long_name = "Precipitation over land",
                    .standard_name = "precipitation_flux",
                    .units = "kg m-2 s-1",
                },
            .size = precipitation_.size(),
        },
        precipitation_);

    coupling_.register_coupled_field(
        registry,
        {
            .role = FieldRole::Export,
            .component = "lnd",
            .attributes =
                {
                    .name = "surface_flux",
                    .long_name = "Surface energy flux",
                    .standard_name = "",
                    .units = "W m-2",
                },
            .size = surface_flux_.size(),
        },
        surface_flux_);
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

} // namespace e3sm::coupler::test

#endif
