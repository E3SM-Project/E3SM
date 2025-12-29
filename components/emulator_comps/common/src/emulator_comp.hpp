/**
 * @file emulator_comp.hpp
 * @brief Abstract base class for all emulated E3SM components.
 *
 * Defines the common interface and shared functionality for emulator
 * components (atmosphere, ocean, ice, land). Derived classes implement
 * component-specific physics through AI/ML inference backends.
 */

#ifndef EMULATOR_COMP_HPP
#define EMULATOR_COMP_HPP

#include "emulator_logger.hpp"
#include <mpi.h>
#include <string>
#include <vector>

namespace emulator {

/**
 * @brief Enumeration of component types in E3SM.
 */
enum class CompType {
  ATM = 0, ///< Atmosphere component
  OCN = 1, ///< Ocean component
  ICE = 2, ///< Sea ice component
  LND = 3  ///< Land component
};

/**
 * @brief Abstract base class for all emulated E3SM components.
 *
 * Provides the common infrastructure for emulator components including:
 * - MPI communicator management
 * - Grid data distribution and storage
 * - Coupling data buffer management
 * - Time stepping and state management
 * - Logging infrastructure
 *
 * Derived classes (e.g., EmulatorAtm) implement the pure virtual methods
 * to provide component-specific behavior.
 *
 * ## Lifecycle
 * 1. `create_instance()` - Initialize MPI and component metadata
 * 2. `set_grid_data()` - Receive grid decomposition from driver
 * 3. `setup_coupling()` - Receive coupling buffer pointers
 * 4. `initialize()` - Load model, read ICs, set up inference
 * 5. `run()` - Execute time step (import → inference → export)
 * 6. `finalize()` - Clean up resources
 *
 * @see EmulatorAtm for the atmosphere implementation
 */
class EmulatorComp {
public:
  /**
   * @brief Construct an emulator component of the given type.
   * @param type Component type (ATM, OCN, ICE, LND)
   */
  explicit EmulatorComp(CompType type);

  virtual ~EmulatorComp() = default;

  // =========================================================================
  // Public Interface (called from Fortran via C interface)
  // =========================================================================

  /**
   * @brief Initialize the component instance.
   *
   * Sets up MPI communicator, stores component metadata, and prepares
   * for grid data and coupling setup.
   *
   * @param comm MPI communicator for this component
   * @param comp_id Component ID assigned by the driver
   * @param input_file Path to component configuration file
   * @param run_type Run type (startup, continue, branch)
   * @param start_ymd Start date as YYYYMMDD
   * @param start_tod Start time of day in seconds
   */
  void create_instance(MPI_Comm comm, int comp_id, const char *input_file,
                       int run_type, int start_ymd, int start_tod);

  /**
   * @brief Set grid decomposition data from the driver.
   *
   * Receives the local portion of the global grid as distributed by
   * the E3SM driver's domain decomposition.
   *
   * @param nx Number of grid points in x-direction (longitude)
   * @param ny Number of grid points in y-direction (latitude)
   * @param num_local_cols Number of columns on this MPI rank
   * @param num_global_cols Total number of columns globally
   * @param col_gids Global IDs for each local column
   * @param lat Latitude values for each local column [radians]
   * @param lon Longitude values for each local column [radians]
   * @param area Grid cell areas for each local column [km²]
   */
  void set_grid_data(int nx, int ny, int num_local_cols, int num_global_cols,
                     const int *col_gids, const double *lat, const double *lon,
                     const double *area);

  /**
   * @brief Set up coupling data buffers.
   *
   * Receives pointers to the MCT attribute vector buffers for data
   * exchange with other components via the coupler.
   *
   * @param import_data Pointer to import buffer (x2a fields)
   * @param export_data Pointer to export buffer (a2x fields)
   * @param num_imports Number of import fields
   * @param num_exports Number of export fields
   * @param field_size Size of each field (should equal num_local_cols)
   */
  void setup_coupling(double *import_data, double *export_data, int num_imports,
                      int num_exports, int field_size);

  /**
   * @brief Initialize the component (phase 2).
   *
   * Called after grid and coupling setup. Loads AI model, reads initial
   * conditions, and prepares for time stepping.
   */
  void initialize();

  /**
   * @brief Execute one time step.
   *
   * Performs the main computational cycle:
   *
   * 1. Import fields from coupler
   * 2. Run AI inference
   * 3. Export fields to coupler
   *
   * @param dt Time step size in seconds
   */
  void run(int dt);

  /**
   * @brief Finalize and clean up the component.
   *
   * Releases resources, finalizes inference backend, and performs
   * any necessary cleanup.
   */
  void finalize();

  // =========================================================================
  // Grid Accessors
  // =========================================================================

  /** @brief Get number of columns on this MPI rank. */
  int get_num_local_cols() const { return m_num_local_cols; }

  /** @brief Get total number of columns globally. */
  int get_num_global_cols() const { return m_num_global_cols; }

  /** @brief Get number of grid points in x-direction (longitude). */
  int get_nx() const { return m_nx; }

  /** @brief Get number of grid points in y-direction (latitude). */
  int get_ny() const { return m_ny; }

  /**
   * @brief Get global IDs for local columns.
   * @param gids Output array (must be pre-allocated with num_local_cols)
   */
  void get_local_col_gids(int *gids) const;

  /**
   * @brief Get lat/lon for local columns.
   * @param lat Output latitude array [radians]
   * @param lon Output longitude array [radians]
   */
  void get_cols_latlon(double *lat, double *lon) const;

  /**
   * @brief Get area for local columns.
   * @param area Output area array [km²]
   */
  void get_cols_area(double *area) const;

  // =========================================================================
  // Component Accessors
  // =========================================================================

  /** @brief Get component type. */
  CompType type() const { return m_type; }

  /** @brief Get MPI communicator. */
  MPI_Comm comm() const { return m_comm; }

  /** @brief Get component ID. */
  int comp_id() const { return m_comp_id; }

  /** @brief Get MPI rank within component communicator. */
  int rank() const { return m_rank; }

  /** @brief Check if this is the root rank (rank 0). */
  bool is_root() const { return m_rank == 0; }

protected:
  // =========================================================================
  // Virtual Methods (implement in derived classes)
  // =========================================================================

  /** @brief Component-specific initialization. */
  virtual void init_impl() = 0;

  /** @brief Component-specific time step execution. */
  virtual void run_impl(int dt) = 0;

  /** @brief Component-specific finalization. */
  virtual void final_impl() = 0;

  /**
   * @brief Run AI inference on packed input/output vectors.
   * @param inputs Packed input features
   * @param outputs Packed output features (will be resized as needed)
   */
  virtual void run_inference(const std::vector<double> &inputs,
                             std::vector<double> &outputs) = 0;

  /** @brief Import fields from coupler (override as needed). */
  virtual void import_from_coupler() {}

  /** @brief Export fields to coupler (override as needed). */
  virtual void export_to_coupler() {}

  // =========================================================================
  // Protected Data
  // =========================================================================

  MPI_Comm m_comm; ///< MPI communicator
  int m_comp_id;   ///< Component ID from driver
  int m_rank;      ///< MPI rank within component
  int m_nprocs;    ///< Number of MPI processes
  CompType m_type; ///< Component type enum
  int m_run_type;  ///< Run type (startup/continue/branch)
  Logger m_logger; ///< Component logger

  // Time and Steps
  int m_current_ymd = 0; ///< Current date as YYYYMMDD
  int m_current_tod = 0; ///< Current time of day [seconds]
  int m_step_count = 0;  ///< Number of steps executed

  // Grid data
  int m_num_local_cols = 0;    ///< Local column count
  int m_num_global_cols = 0;   ///< Global column count
  int m_nx = 0;                ///< Grid points in x (longitude)
  int m_ny = 0;                ///< Grid points in y (latitude)
  std::vector<int> m_col_gids; ///< Global IDs for local columns
  std::vector<double> m_lat;   ///< Latitude [radians]
  std::vector<double> m_lon;   ///< Longitude [radians]
  std::vector<double> m_area;  ///< Cell area [km²]

  // Coupling data
  double *m_import_data = nullptr; ///< Import buffer pointer (x2a)
  double *m_export_data = nullptr; ///< Export buffer pointer (a2x)
  int m_num_imports = 0;           ///< Number of import fields
  int m_num_exports = 0;           ///< Number of export fields
  int m_field_size = 0;            ///< Size per field (should = num_local_cols)

  std::string m_input_file; ///< Path to configuration file

private:
  void read_grid_file(const std::string &config_file);
  void distribute_grid_data(const std::vector<double> &lon_global,
                            const std::vector<double> &lat_global,
                            const std::vector<double> &area_global);
  void setup_default_grid();
  void advance_time(int dt);

  bool m_initialized = false;
};

// ===========================================================================
// Utility Functions
// ===========================================================================

/**
 * @brief Get the coupling import prefix for a component type.
 * @param type Component type
 * @return Prefix string (e.g., "x2a" for atmosphere)
 */
inline std::string get_import_prefix(CompType type) {
  switch (type) {
  case CompType::ATM:
    return "x2a";
  case CompType::OCN:
    return "x2o";
  case CompType::ICE:
    return "x2i";
  case CompType::LND:
    return "x2l";
  default:
    return "x2x";
  }
}

/**
 * @brief Get the coupling export prefix for a component type.
 * @param type Component type
 * @return Prefix string (e.g., "a2x" for atmosphere)
 */
inline std::string get_export_prefix(CompType type) {
  switch (type) {
  case CompType::ATM:
    return "a2x";
  case CompType::OCN:
    return "o2x";
  case CompType::ICE:
    return "i2x";
  case CompType::LND:
    return "l2x";
  default:
    return "x2x";
  }
}

/**
 * @brief Get the short name for a component type.
 * @param type Component type
 * @return Short name string (e.g., "atm" for atmosphere)
 */
inline std::string get_comp_name(CompType type) {
  switch (type) {
  case CompType::ATM:
    return "atm";
  case CompType::OCN:
    return "ocn";
  case CompType::ICE:
    return "ice";
  case CompType::LND:
    return "lnd";
  default:
    return "xxx";
  }
}

} // namespace emulator

#endif // EMULATOR_COMP_HPP
