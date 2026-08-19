#ifndef E3SM_COUPLER_TESTS_FAKE_COMPONENTS_HPP
#define E3SM_COUPLER_TESTS_FAKE_COMPONENTS_HPP

#include <algorithm>
#include <field_registry.hpp>
#include <string_view>
#include <vector>

namespace e3sm::coupler::test {

struct FieldBinding {
  FieldID id;
  std::span<double> data;
};

class CouplingAdapter {
public:
  void register_coupled_field(FieldRegistry& registry, RegisteredField metadata,
                              std::span<double> data) {
    const auto id = registry.register_field(std::move(metadata));

    fields_.push_back({
        .id = id,
        .data = data,
    });
  }

  void configure(std::span<const FieldID> export_ids,
                 std::span<const FieldID> import_ids) {
    active_exports_ = resolve(export_ids);
    active_imports_ = resolve(import_ids);
  }

  void export_fields(std::span<ExportBuffer> buffers) const {
    if (buffers.size() != active_exports_.size()) {
      throw std::runtime_error("Incorrect number of export buffers");
    }

    for (std::size_t i = 0; i < buffers.size(); ++i) {
      const auto& local = active_exports_[i];
      auto& remote = buffers[i];

      if (local.id != remote.id) {
        throw std::runtime_error("Export field ID mismatch");
      }

      if (local.data.size() != remote.data.size()) {
        throw std::runtime_error("Export field size mismatch");
      }

      std::ranges::copy(local.data, remote.data.begin());
    }
  }

  void import_fields(std::span<const ImportBuffer> buffers) {
    if (buffers.size() != active_imports_.size()) {
      throw std::runtime_error("Incorrect number of import buffers");
    }

    for (std::size_t i = 0; i < buffers.size(); ++i) {
      auto& local = active_imports_[i];
      const auto& remote = buffers[i];

      if (local.id != remote.id) {
        throw std::runtime_error("Import field ID mismatch");
      }

      if (local.data.size() != remote.data.size()) {
        throw std::runtime_error("Import field size mismatch");
      }

      std::ranges::copy(remote.data, local.data.begin());
    }
  }

private:
  /**
   * @brief Helper to resolve the FieldBinding corresponding to a FieldID
   * */
  std::vector<FieldBinding> resolve(std::span<const FieldID> ids) const {
    std::vector<FieldBinding> result;
    result.reserve(ids.size());

    for (FieldID id : ids) {
      auto it = std::ranges::find_if(
          fields_, [id](const FieldBinding& field) { return field.id == id; });

      if (it == fields_.end()) {
        throw std::runtime_error(
            "Coupling configuration references unknown field");
      }

      result.push_back(*it);
    }

    return result;
  }

  std::vector<FieldBinding> fields_;
  std::vector<FieldBinding> active_exports_;
  std::vector<FieldBinding> active_imports_;
};

// FAKE COMPONENTS

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

  void configure_coupling(std::span<const FieldID> export_ids,
                          std::span<const FieldID> import_ids) {
    coupling_.configure(export_ids, import_ids);
  }

  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }

  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }

  void run() {
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

  void run() {
    for (std::size_t i = 0; i < surface_flux_.size(); ++i) {
      surface_flux_[i] = 2.0 * precipitation_[i];
    }
  }

  void configure_coupling(std::span<const FieldID> export_ids,
                          std::span<const FieldID> import_ids) {
    coupling_.configure(export_ids, import_ids);
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
