#ifndef E3SM_COUPLER_COMPONENTS_HPP
#define E3SM_COUPLER_COMPONENTS_HPP
#include <algorithm>
#include <concepts>
#include <field_registry.hpp>
#include <span>
#include <string>
#include <string_view>

namespace e3sm::coupler {

template <typename T>
concept E3SMComponent =
    requires(T& component, std::string_view name, FieldRegistry& registry,
             std::span<FieldID> export_ids, std::span<FieldID> import_ids,
             std::span<ExportBuffer> export_buffer,
             std::span<const ImportBuffer> import_buffer) {
      { component.name() } -> std::convertible_to<std::string_view>;
      { component.populate_registry(registry) } -> std::same_as<void>;
      {
        component.configure_coupling(export_ids, import_ids)
      } -> std::same_as<void>;
      { component.run() } -> std::same_as<void>;
      { component.export_fields(export_buffer) } -> std::same_as<void>;
      { component.import_fields(import_buffer) } -> std::same_as<void>;
    };

/**
 * @brief Class to hold any object that satifies `E3SMComponent` to be used by
 * the CouplerDriver
 */
class AnyComponent {
public:
  template <E3SMComponent C>
  explicit AnyComponent(C& component)
      : object_(&component),
        name_([](void* obj) { return static_cast<C*>(obj)->name(); }),
        populate_registry_([](void* obj, FieldRegistry& registry) {
          static_cast<C*>(obj)->populate_registry(registry);
        }),
        configure_coupling_([](void* obj, std::span<FieldID> export_ids,
                               std::span<FieldID> import_ids) {
          static_cast<C*>(obj)->configure_coupling(export_ids, import_ids);
        }),
        run_([](void* obj) { static_cast<C*>(obj)->run(); }),
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
  void configure_coupling(std::span<FieldID> export_ids,
                          std::span<FieldID> import_ids) {
    configure_coupling_(object_, export_ids, import_ids);
  }
  void run() { run_(object_); }
  void export_fields(std::span<ExportBuffer> buf) {
    export_fields_(object_, buf);
  }
  void import_fields(std::span<ImportBuffer> buf) {
    import_fields_(object_, buf);
  }

private:
  using NameFn = std::string_view (*)(void*);
  using PopulateRegistryFn = void (*)(void*, FieldRegistry&);
  using ConfigCouplingFn = void (*)(void*, std::span<FieldID>,
                                    std::span<FieldID>);
  using RunFn = void (*)(void*);
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
} // namespace e3sm::coupler
#endif
