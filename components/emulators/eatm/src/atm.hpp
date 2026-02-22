/**
 * @file atm.hpp
 * @brief Atmosphere emulator component declaration.
 *
 * Defines the EmulatorAtm class which implements an AI-based atmosphere
 * component for E3SM. Inherits from the Emulator base class and adds
 * atmosphere-specific coupling, field management, and inference.
 */

#ifndef EMULATOR_ATM_HPP
#define EMULATOR_ATM_HPP

#include "emulator.hpp"
#include <memory>
#include <string>
#include <vector>

namespace emulator {

/**
 * @brief Atmosphere emulator component.
 *
 * Derived from Emulator, provides atmosphere-specific functionality:
 * - Coupling field mappings (x2a inputs, a2x outputs)
 * - AI model integration via configurable inference backends
 * - MCT interface for CIME integration
 *
 * Currently assumes a structured lat-lon grid. Grid dimensions (nx, ny)
 * are read from atm_in and used to compute the total global column count
 * as nx * ny. Lat/lon coordinates are stored and passed to MCT in degrees.
 *
 * ## Lifecycle
 * 1. Constructor creates EmulatorAtm with ATM_COMP type
 * 2. create_instance() sets MPI, comp_id, parses config for grid dims
 * 3. set_grid_data() sets spatial decomposition (optional override)
 * 4. init_coupling_indices() parses MCT field lists
 * 5. setup_coupling() sets buffer pointers
 * 6. initialize() loads model and reads initial conditions
 * 7. run() executes time steps (import -> inference -> export)
 * 8. finalize() cleans up resources
 */
class EmulatorAtm : public Emulator {
public:
  EmulatorAtm();
  ~EmulatorAtm() override = default;

  // =========================================================================
  // Setup methods (called before initialize)
  // =========================================================================

  /**
   * @brief Set MPI communicator, component ID, and run settings.
   */
  void create_instance(int comm, int comp_id,
                       const std::string &input_file,
                       int run_type, int start_ymd, int start_tod);

  /**
   * @brief Set grid decomposition data from driver.
   */
  void set_grid_data(int nx, int ny,
                     int num_local_cols, int num_global_cols,
                     const int *col_gids,
                     const double *lat, const double *lon,
                     const double *area);

  /**
   * @brief Initialize coupling field indices from MCT field lists.
   */
  void init_coupling_indices(const std::string &export_fields,
                             const std::string &import_fields);

  /**
   * @brief Set up coupling buffer pointers from MCT.
   */
  void setup_coupling(double *import_data, double *export_data,
                      int num_imports, int num_exports,
                      int field_size);

  // =========================================================================
  // Accessors
  // =========================================================================

  int get_num_local_cols() const { return m_num_local_cols; }
  int get_num_global_cols() const { return m_num_global_cols; }
  int get_nx() const { return m_nx; }
  int get_ny() const { return m_ny; }
  void get_local_col_gids(int *gids) const;
  void get_cols_latlon(double *lat, double *lon) const;
  void get_cols_area(double *area) const;

protected:
  // Virtual methods from Emulator base
  void init_impl() override;
  void run_impl(int dt) override;
  void final_impl() override;

private:
  // =========================================================================
  // Grid and decomposition
  // =========================================================================
  int m_nx = 0;                ///< Grid x-dimension
  int m_ny = 0;                ///< Grid y-dimension
  int m_num_local_cols = 0;    ///< Local columns on this rank
  int m_num_global_cols = 0;   ///< Total global columns
  std::vector<int> m_col_gids; ///< Global IDs for local columns
  std::vector<double> m_lat;   ///< Latitude [degrees]
  std::vector<double> m_lon;   ///< Longitude [degrees]
  std::vector<double> m_area;  ///< Cell areas

  // =========================================================================
  // Coupling
  // =========================================================================
  double *m_import_data = nullptr; ///< MCT import buffer pointer
  double *m_export_data = nullptr; ///< MCT export buffer pointer
  int m_num_imports = 0;           ///< Number of import fields
  int m_num_exports = 0;           ///< Number of export fields

  // =========================================================================
  // Configuration
  // =========================================================================
  int m_comm = 0;              ///< MPI communicator
  std::string m_input_file;    ///< Path to atm_in config file
  int m_run_type = 0;          ///< Run type (startup/continue/branch)

  // =========================================================================
  // Helper methods
  // =========================================================================
  void import_coupling_fields();
  void export_coupling_fields();
  void prepare_inputs();
  void process_outputs();
};

} // namespace emulator

#endif // EMULATOR_ATM_HPP
