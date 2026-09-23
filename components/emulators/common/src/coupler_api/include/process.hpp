#ifndef E3SM_COUPLER_API_PROCESS_HPP
#define E3SM_COUPLER_API_PROCESS_HPP
#include "coupler.hpp"
#include <algorithm>
#include <field_registry.hpp>
#include <span>
#include <string_view>
#include <vector>

namespace e3sm::coupler {

struct FieldBinding {
  FieldID id;
  std::span<double> data;
};

struct ProcessRunOpts {
  float dt;
};

template <typename T>
concept CoupledProcess =
    requires(T& component, std::string_view name, FieldRegistry& registry,
             const ActiveCouplingFields& plan, const ProcessRunOpts run_opts,
             std::span<ExportBuffer> export_buffer,
             std::span<const ImportBuffer> import_buffer) {
      { component.name() } -> std::convertible_to<std::string_view>;
      { component.populate_registry(registry) } -> std::same_as<void>;
      { component.configure_coupling(plan) } -> std::same_as<void>;
      { component.run(run_opts) } -> std::same_as<void>;
      { component.export_fields(export_buffer) } -> std::same_as<void>;
      { component.import_fields(import_buffer) } -> std::same_as<void>;
    };

/**
 * @brief Class to hold any object that satifies `CoupledProcess` to be used by
 * the CouplerDriver
 */
class AnyComponent {
public:
  template <CoupledProcess C>
  explicit AnyComponent(C& component)
      : object_(&component),
        name_([](void* obj) { return static_cast<C*>(obj)->name(); }),
        populate_registry_([](void* obj, FieldRegistry& registry) {
          static_cast<C*>(obj)->populate_registry(registry);
        }),
        configure_coupling_([](void* obj, const ActiveCouplingFields& plan) {
          static_cast<C*>(obj)->configure_coupling(plan);
        }),
        run_([](void* obj, const ProcessRunOpts run_opts) { static_cast<C*>(obj)->run(run_opts); }),
        export_fields_([](void* obj, std::span<ExportBuffer> buf) {
          static_cast<C*>(obj)->export_fields(buf);
        }),
        import_fields_([](void* obj, std::span<ImportBuffer> buf) {
          static_cast<C*>(obj)->import_fields(buf);
        }) {}

  std::string_view name() const { return name_(object_); }
  void populate_registry(FieldRegistry& registry) {
    populate_registry_(object_, registry);
  }
  void configure_coupling(const ActiveCouplingFields& plan) {
    configure_coupling_(object_, plan);
  }
  void run(const ProcessRunOpts run_opts) { run_(object_, run_opts); }
  void export_fields(std::span<ExportBuffer> buf) {
    export_fields_(object_, buf);
  }
  void import_fields(std::span<ImportBuffer> buf) {
    import_fields_(object_, buf);
  }

private:
  using NameFn = std::string_view (*)(void*);
  using PopulateRegistryFn = void (*)(void*, FieldRegistry&);
  using ConfigCouplingFn = void (*)(void*, const ActiveCouplingFields&);
  using RunFn = void (*)(void*, const ProcessRunOpts);
  using ExportFn = void (*)(void*, std::span<ExportBuffer>);
  using ImportFn = void (*)(void*, std::span<ImportBuffer>);

  void* object_;
  NameFn name_;
  PopulateRegistryFn populate_registry_;
  ConfigCouplingFn configure_coupling_;
  RunFn run_;
  ExportFn export_fields_;
  ImportFn import_fields_;
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

  void configure(const ActiveCouplingFields& plan) {
    active_exports_ = resolve(plan.export_ids);
    active_imports_ = resolve(plan.import_ids);
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
} // namespace e3sm::coupler
#endif
