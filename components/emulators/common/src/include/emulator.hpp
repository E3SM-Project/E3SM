/**
 * @file emulator.hpp
 * @brief Class for E3SM emulators
 */
#ifndef E3SM_EMULATOR_HPP
#define E3SM_EMULATOR_HPP

#include "moab_adapter.hpp"
#include "moab_transport.hpp"
#include <coupler.hpp>
#include <emulator_config.hpp>
#include <field_registry.hpp>
#include <inference_backend.hpp>
#include <memory>
#include <mpi.h>
#include <optional>
#include <process.hpp>

#include "grid_types.hpp"
#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace e3sm::emulator {

using namespace e3sm::coupler;

/**
 * @brief Abstract base class for all E3SM emulators.
 *
 * Provides the common infrastructure for emulators.
 * Derived classes implement the pure virtual methods for
 * emulator-specific behavior.
 */
class Emulator {
public:
  Emulator(const GridDesc grid, const EmulatorConfig emu_config);
  /// The backend and MOAB adapter clean up after themselves. Call finalize()
  /// to see a failure in the model's finalize(), which the destructor can only
  /// log. With MOAB, destroy the emulator before iMOAB_Finalize and MPI_Finalize.
  ~Emulator() = default;

  // Coupled spans point into this object's field buffers
  Emulator(const Emulator&) = delete;
  Emulator& operator=(const Emulator&) = delete;
  Emulator(Emulator&&) = delete;
  Emulator& operator=(Emulator&&) = delete;

  // Lifecycle methods
  void initialize();
  void run(const ProcessRunOpts opts);
  void finalize();

  // Accessors
  std::string_view name() const { return m_name; }
  MPI_Comm comm() const { return m_grid.comm; }
  EmulatorType type() const { return m_type; }
  int id() const { return m_component_id; }
  int moab_app_id() const { return m_moab_app_id; }
  bool is_initialized() const { return m_initialized; }
  int step_count() const { return m_step_count; }
  void print_info(std::ostream& os) const;

  // Common coupling interface
  void populate_registry(FieldRegistry& registry);

  void configure_coupling(const ActiveCouplingFields& plan) {
    m_coupling.configure(plan);
  }

  void export_fields(std::span<ExportBuffer> buffers) const {
    m_coupling.export_fields(buffers);
  }

  void import_fields(std::span<const ImportBuffer> buffers) {
    m_coupling.import_fields(buffers);
  }
  std::span<const int> get_local_gids() const {
    return m_coupling.get_local_gids();
  };

  // TODO: Let's move this to protected after adjusting Fortran-C API??
  void set_grid_data(const GridDesc& grid);
  int get_num_local_cols() const { return m_grid.num_local_columns; };
  int get_num_global_cols() const { return m_grid.num_global_columns; };
  int get_nx() const { return m_grid.nx; };
  int get_ny() const { return m_grid.ny; };
  std::span<const double> get_latitudes() const { return m_grid.lat; };
  std::span<const double> get_longitudes() const { return m_grid.lon; };
  std::span<const double> get_areas() const { return m_grid.area; };

  // diagnostics and lookups
  std::span<const EmulatorFieldDesc> field_descriptors() const;
  // These throw if the name is both an input and an output; pass a role instead.
  const EmulatorFieldDesc& get_field_descriptor(std::string_view name) const;
  std::span<const double> get_field(std::string_view name);

  const EmulatorFieldDesc& get_field_descriptor(std::string_view name,
                                                ModelFieldRole role) const;
  std::span<const double> get_field(std::string_view name, ModelFieldRole role);
  std::span<double> input_fields() { return std::span<double>{input_fields_}; };
  std::span<double> output_fields() {
    return std::span<double>{output_fields_};
  };

protected:
  template <typename Container>
  void register_import_field(FieldRegistry& registry, RegisteredField field,
                             Container& storage) {
    field.role = FieldRole::Import;
    field.component = m_name;
    field.size = storage.size();

    m_coupling.register_coupled_field(registry, std::move(field), storage);
  }
  template <typename Container>
  void register_export_field(FieldRegistry& registry, RegisteredField field,
                             Container& storage) {
    field.role = FieldRole::Export;
    field.component = m_name;
    field.size = storage.size();

    m_coupling.register_coupled_field(registry, std::move(field), storage);
  }

  // Returns a view to part of the flat buffer that holds the field's data
  std::span<double> field_view(const EmulatorFieldDesc& field);
  std::span<double> field_view(std::string_view name);

  // Call this to register the emulator as an application with MOAB
  void register_with_moab();

private:
  std::filesystem::path m_config_path;
  GridDesc m_grid;
  EmulatorConfig m_emulator_config;
  CouplingAdapter m_coupling;
  EmulatorType m_type;
  int m_component_id;
  std::string m_name;
  int m_moab_app_id;
  std::optional<MoabAdapter> m_moab;

  bool m_initialized = false;
  bool m_finalized = false;
  int m_step_count = 0;

  // dimensions_ holds the domain info over which the emulator needs to be
  // allocated.
  const DimensionDesc& dimension(std::string_view name) const;
  // initalization helpers
  void apply_config(const EmulatorConfig& config);

  std::vector<EmulatorFieldDesc> field_descriptors_;
  std::size_t local_field_size(const FieldConfig& field) const;
  RegisteredField make_registered_field(const EmulatorFieldDesc& field);
  std::size_t append_fields(std::span<const std::string> order,
                            std::span<const FieldConfig>, ModelFieldRole kind);

  MoabTagDescriptor make_moab_tag(const EmulatorFieldDesc& desc) const;
  // Flat Buffers for inputs and outputs
  std::vector<double> input_fields_;
  std::vector<double> output_fields_;

  std::unique_ptr<inference::InferenceBackend> backend_;
};

} // namespace e3sm::emulator

#endif // E3SM_EMULATOR_HPP
