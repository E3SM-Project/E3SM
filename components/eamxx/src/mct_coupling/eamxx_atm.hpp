/**
 * @file atm.hpp
 * @brief Atmosphere emulator component declaration.
 *
 * Defines the Atm class which implements an AI-based atmosphere
 * component for E3SM. Inherits from the Emulator base class and adds
 * atmosphere-specific coupling, field management, and inference.
 */

#ifndef EAMXX_ATM_HPP
#define EAMXX_ATM_HPP

#include <emulator.hpp> // should be called "component"
#include <emulator_c_api.hpp>
#include <memory>
#include <string>
#include <vector>

#include <mpi.h>

namespace eamxx {

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
 * 1. Constructor creates Atm with ATM_COMP type
 * 2. create_instance() sets MPI, comp_id, parses config for grid dims
 * 3. set_grid_data() sets spatial decomposition (optional override)
 * 4. init_coupling_indices() parses MCT field lists
 * 5. setup_coupling() sets buffer pointers
 * 6. initialize() loads model and reads initial conditions
 * 7. run() executes time steps (import -> inference -> export)
 * 8. finalize() cleans up resources
 */
class Atm : public Emulator {
public:
  Atm();

  Atm(const Atm&) = delete;
  ~Atm() = default;

  Atm& operator=(const Atm&) = delete;

  // =========================================================================
  // Setup methods (called before initialize)
  // =========================================================================

  /**
   * @brief Create and register an EAMxx instance with the given options.
   */
  void create_instance(const MPI_Fint f_comm, const int atm_id,
                       const char* input_yaml_file,
                       const char* atm_log_file,
                       const int run_type,
                       const int run_start_ymd,
                       const int run_start_tod,
                       const int case_start_ymd,
                       const int case_start_tod,
                       const char* calendar_name);{
                     /*const char* caseid,
                       const char* rest_caseid,
                       const char* hostname,
                       const char* username,
                       const char* versionid);*/
  /**
   * @brief Set grid decomposition data from driver.
   */
  // FIXME: does this pertain to EAMxx?
  void set_grid_data(const EmulatorGridDesc& grid) override;

  /**
   * @brief Initialize coupling field indices from MCT field lists.
   */
  void init_coupling_indices(const std::string &export_fields,
                             const std::string &import_fields) override;

  /**
   * @brief Set up coupling buffer pointers from MCT.
   */
  void setup_coupling(const EmulatorCouplingDesc& cpl) override;

  // =========================================================================
  // Accessors
  // =========================================================================

  int get_num_local_cols() const override;
  int get_num_global_cols() const override;
  void get_local_col_gids(int *gids) const override;
  void get_cols_latlon(double *lat, double *lon) const override;
  void get_cols_area(double *area) const override;

protected:
  // Virtual methods from Emulator base
  void init_impl() override;
  void run_impl(int dt) override;
  void final_impl() override;
  void print_extra_info(std::ostream& os) const override {};

private:
  
  // =========================================================================
  // Grid and decomposition
  // =========================================================================
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
  std::string m_log_file;      ///< Path to log file
  int m_run_type = 0;          ///< Run type (startup/continue/branch)

  // =========================================================================
  // Helper methods
  // =========================================================================
  void import_coupling_fields();
  void export_coupling_fields();
  void prepare_inputs();
  void process_outputs();
};

} // namespace eamxx

#endif // EAMXX_ATM_HPP
