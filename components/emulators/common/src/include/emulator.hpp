/**
 * @file emulator.hpp
 * @brief Abstract base class for all E3SM emulators.
 */
#ifndef E3SM_EMULATOR_HPP
#define E3SM_EMULATOR_HPP

#include <cstddef>
#include <mpi.h>

#include "coupler_types.hpp"
#include "emulator_c_api.hpp"
#include <coupler.hpp>
#include <emulator_config.hpp>
#include <field_registry.hpp>
#include <process.hpp>
#include <span>
#include <string>
#include <vector>

namespace emulator {

using namespace e3sm::coupler;

/**
 * @brief: Struct that describes how fields is stored in emulator and its
 purpose
 * @fields:
    - ModelFieldKind role;
    - FieldConfig config;
    - std::size_t offset;
    - std::size_t size;
 */
struct EmulatorFieldDesc {
  ModelFieldRole role;
  FieldConfig config;
  std::size_t offset;
  std::size_t size;
};

/**
 * @brief C++-version EmulatorDimensionDesc
  -std::string name;
  -std::size_t local_extent;
  -std::size_t global_extent;
 */
struct DimensionDesc {
  std::string name;
  std::size_t local_extent;
  std::size_t global_extent;
};
std::string to_string(const DimensionDesc& desc);

/**
 * @brief Struct for holding the grid info for this process
 *  and it's possible dimentions;
   @fields:
   - int grid_type;
   - int nx;
   - int ny;
   - std::vector<DimensionDesc> dimensions;

   - std::vector<int> col_gids;
   - std::vector<double> lat;
   - std::vector<double> lon;
   - std::vector<double> area;
 *
 */
struct GridDesc {
  int grid_type;
  int nx;
  int ny;
  int num_global_columns;
  int num_local_columns;
  MPI_Comm comm;

  std::vector<DimensionDesc> dimensions;
  std::vector<int> col_gids;
  std::vector<double> lat;
  std::vector<double> lon;
  std::vector<double> area;
};

/**
 * @brief:
  std::filesystem::path path;
  GridDesc grid;
 */
struct CreateConfig {
  std::filesystem::path path;
  GridDesc grid;
};

/**
 * @brief Abstract base class for all E3SM emulators.
 *
 * Provides the common infrastructure for emulators.
 * Derived classes implement the pure virtual methods for
 * emulator-specific behavior.
 */
class Emulator {
public:
  /**
   * @brief Construct a new Emulator.
   *
   * @param comm MPI communicator
   * @param id Emulator ID (-1 if unassigned)
   */
  Emulator(const CreateConfig);
  virtual ~Emulator() = default;

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

  // TODO: Let's move this to protected after adjusting Fortran-C API??
  void set_grid_data(const EmulatorGridDesc& grid);
  int get_num_local_cols() const { return m_grid.num_local_columns; };
  int get_num_global_cols() const { return m_grid.num_global_columns; };
  int get_nx() const { return m_grid.nx; };
  int get_ny() const { return m_grid.ny; };
  std::span<const int> get_local_col_gids() const { return m_grid.col_gids; };
  std::span<const double> get_latitudes() const { return m_grid.lat; };
  std::span<const double> get_longitudes() const { return m_grid.lon; };
  std::span<const double> get_areas() const { return m_grid.area; };

  // diagnostics and lookups
  std::span<const EmulatorFieldDesc> field_descriptors() const;
  const EmulatorFieldDesc& get_field_descriptor(std::string_view name) const;
  std::span<const double> get_field(std::string_view name);
  std::span<double> input_fields() { return std::span<double>{input_fields_}; };
  std::span<double> output_fields() {
    return std::span<double>{output_fields_};
  };

  // Emulator specific field initialization
  virtual void init_data() = 0;
  // Emulator specific run details
  virtual void run_impl(const ProcessRunOpts opts) = 0;
  // Emulator specific diagnostics
  virtual void print_extra_info(std::ostream& os) const {}

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
  CouplingAdapter m_coupling;
  EmulatorType m_type;
  int m_component_id;
  std::string m_name;
  int m_moab_app_id;

  bool m_initialized = false;
  int m_step_count = 0;

  // dimensions_ holds the domain info over which the emulator needs to be
  // allocated.
  const DimensionDesc& dimension(std::string_view name) const;

  std::vector<EmulatorFieldDesc> field_descriptors_;
  std::size_t local_field_size(const FieldConfig& field) const;
  RegisteredField make_registered_field(const EmulatorFieldDesc& field);
  std::size_t append_fields(std::span<const FieldConfig>, ModelFieldRole kind);

  // Flat Buffers for inputs and outputs
  std::vector<double> input_fields_;
  std::vector<double> output_fields_;
};

} // namespace emulator

#endif // E3SM_EMULATOR_HPP
