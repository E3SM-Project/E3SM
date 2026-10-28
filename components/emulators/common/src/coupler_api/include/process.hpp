#ifndef E3SM_COUPLER_API_PROCESS_HPP
#define E3SM_COUPLER_API_PROCESS_HPP
#include "coupler.hpp"
#include <field_registry.hpp>
#include <memory>
#include <span>
#include <string_view>
#include <typeindex>
#include <vector>

namespace e3sm::coupler {

struct FieldBinding {
  FieldID id;
  std::span<double> data;
};

struct ProcessRunOpts {
  float dt;
  int coupling_interval = 1;
};

template <typename T>
concept CoupledProcess =
    requires(T& component, std::string_view name, FieldRegistry& registry,
             const ActiveCouplingFields& plan, const ProcessRunOpts run_opts,
             std::span<ExportBuffer> export_buffers,
             std::span<const ImportBuffer> import_buffers) {
      { component.name() } -> std::convertible_to<std::string_view>;
      { component.populate_registry(registry) } -> std::same_as<void>;
      { component.configure_coupling(plan) } -> std::same_as<void>;
      { component.run(run_opts) } -> std::same_as<void>;
      { component.export_fields(export_buffers) } -> std::same_as<void>;
      { component.import_fields(import_buffers) } -> std::same_as<void>;
    };

/**
 * @brief Class to hold any object that satifies `CoupledProcess` to be used by
 * the CouplerDriver
 */
class AnyComponent {
public:
  template <CoupledProcess C>
  AnyComponent(std::unique_ptr<C> component)
      : object_(component.release(),
                [](void* obj) { delete static_cast<C*>(obj); }),
        type_(typeid(C)),
        name_([](void* obj) { return static_cast<C*>(obj)->name(); }),
        populate_registry_([](void* obj, FieldRegistry& registry) {
          static_cast<C*>(obj)->populate_registry(registry);
        }),
        configure_coupling_([](void* obj, const ActiveCouplingFields& plan) {
          static_cast<C*>(obj)->configure_coupling(plan);
        }),
        run_([](void* obj, const ProcessRunOpts run_opts) {
          static_cast<C*>(obj)->run(run_opts);
        }),
        export_fields_([](void* obj, std::span<ExportBuffer> buf) {
          static_cast<C*>(obj)->export_fields(buf);
        }),
        import_fields_([](void* obj, std::span<ImportBuffer> buf) {
          static_cast<C*>(obj)->import_fields(buf);
        }) {}

  std::string_view name() const { return name_(object_.get()); }
  void populate_registry(FieldRegistry& registry) {
    populate_registry_(object_.get(), registry);
  }
  void configure_coupling(const ActiveCouplingFields& plan) {
    configure_coupling_(object_.get(), plan);
  }
  void run(const ProcessRunOpts run_opts) { run_(object_.get(), run_opts); }
  void export_fields(std::span<ExportBuffer> buf) {
    export_fields_(object_.get(), buf);
  }
  void import_fields(std::span<ImportBuffer> buf) {
    import_fields_(object_.get(), buf);
  }

  template <typename C> C& as() {
    if (type_ != std::type_index(typeid(C))) {
      throw std::bad_cast{};
    }
    return *static_cast<C*>(object_.get());
  }

private:
  using NameFn = std::string_view (*)(void*);
  using PopulateRegistryFn = void (*)(void*, FieldRegistry&);
  using ConfigCouplingFn = void (*)(void*, const ActiveCouplingFields&);
  using RunFn = void (*)(void*, const ProcessRunOpts);
  using ExportFn = void (*)(void*, std::span<ExportBuffer>);
  using ImportFn = void (*)(void*, std::span<ImportBuffer>);

  using DeleteFn = void (*)(void*);
  std::unique_ptr<void, DeleteFn> object_;
  std::type_index type_;
  NameFn name_;
  PopulateRegistryFn populate_registry_;
  ConfigCouplingFn configure_coupling_;
  RunFn run_;
  ExportFn export_fields_;
  ImportFn import_fields_;
};

// helper class for implementing processes / components
class CouplingAdapter {
public:
  void register_coupled_field(FieldRegistry& registry, RegisteredField metadata,
                              std::span<double> data);

  void configure(const ActiveCouplingFields& plan);

  void export_fields(std::span<ExportBuffer> buffers) const;

  void import_fields(std::span<const ImportBuffer> buffers);

  std::vector<FieldBinding> resolve(std::span<const FieldID> ids) const;

  std::vector<FieldBinding> fields_;
  std::vector<FieldBinding> active_exports_;
  std::vector<FieldBinding> active_imports_;
};
} // namespace e3sm::coupler
#endif
