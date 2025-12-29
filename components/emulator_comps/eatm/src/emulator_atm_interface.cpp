/**
 * @file emulator_atm_interface.cpp
 * @brief C interface for the atmosphere emulator (Fortran-callable).
 *
 * Provides extern "C" functions that can be called from Fortran via
 * ISO_C_BINDING. These functions wrap the EmulatorAtm C++ class methods.
 *
 * The interface follows the E3SM component lifecycle:
 * 1. create_instance - Initialize MPI and component
 * 2. set_grid_data - Set grid decomposition
 * 3. init_coupling_indices - Parse MCT field lists
 * 4. setup_coupling - Set coupling buffer pointers
 * 5. init - Load model and prepare for time stepping
 * 6. run - Execute time steps
 * 7. finalize - Clean up resources
 *
 * @see EmulatorAtm for the C++ implementation
 * @see emulator_atm_f2c.F90 for the Fortran interface module
 */

#include "../../common/src/emulator_comp.hpp"
#include "../../common/src/emulator_context.hpp"
#include "emulator_atm.hpp"

namespace {

/**
 * @brief Get non-const reference to the atmosphere emulator.
 * @return Reference to EmulatorAtm instance from context
 */
emulator::EmulatorAtm &get_atm_emulator_nonconst() {
  auto &ctx = emulator::EmulatorContext::singleton();
  return ctx.getNonConst<emulator::EmulatorAtm>();
}

/**
 * @brief Get const reference to the atmosphere emulator.
 * @return Const reference to EmulatorAtm instance from context
 */
const emulator::EmulatorAtm &get_atm_emulator() {
  const auto &ctx = emulator::EmulatorContext::singleton();
  return ctx.get<emulator::EmulatorAtm>();
}

} // anonymous namespace

extern "C" {

// ===========================================================================
// Atmosphere Emulator C Interface (callable from Fortran)
// ===========================================================================

/**
 * @brief Create and initialize the atmosphere emulator instance.
 *
 * @param f_comm Fortran MPI communicator
 * @param comp_id Component ID from driver
 * @param input_file Path to atm_in configuration file
 * @param log_file Path to log file (NULL for stdout)
 * @param run_type Run type (startup, continue, branch)
 * @param start_ymd Start date as YYYYMMDD
 * @param start_tod Start time of day in seconds
 */
void emulator_atm_create_instance(const MPI_Fint f_comm, const int comp_id,
                                  const char *input_file, const char *log_file,
                                  const int run_type, const int start_ymd,
                                  const int start_tod) {
  MPI_Comm c_comm = MPI_Comm_f2c(f_comm);

  auto &ctx = emulator::EmulatorContext::singleton();
  auto &atm = ctx.create<emulator::EmulatorAtm>();

  atm.create_instance(c_comm, comp_id, input_file, run_type, start_ymd,
                      start_tod);
  atm.set_log_file(log_file ? log_file : "");
}

/**
 * @brief Set grid decomposition data from driver.
 */
void emulator_atm_set_grid_data(const int nx, const int ny,
                                const int num_local_cols,
                                const int num_global_cols, const int *col_gids,
                                const double *lat, const double *lon,
                                const double *area) {
  auto &atm = get_atm_emulator_nonconst();
  atm.set_grid_data(nx, ny, num_local_cols, num_global_cols, col_gids, lat, lon,
                    area);
}

/**
 * @brief Initialize coupling field index mappings.
 */
void emulator_atm_init_coupling_indices(const char *export_fields,
                                        const char *import_fields) {
  auto &atm = get_atm_emulator_nonconst();
  std::string exp_str = export_fields ? export_fields : "";
  std::string imp_str = import_fields ? import_fields : "";
  atm.init_coupling_indices(exp_str, imp_str);
}

/**
 * @brief Set up coupling buffer pointers from MCT.
 */
void emulator_atm_setup_coupling(double *import_data, double *export_data,
                                 const int num_imports, const int num_exports,
                                 const int field_size) {
  auto &atm = get_atm_emulator_nonconst();
  atm.setup_coupling(import_data, export_data, num_imports, num_exports,
                     field_size);
}

/**
 * @brief Initialize the atmosphere emulator (phase 2).
 */
void emulator_atm_init() {
  auto &atm = get_atm_emulator_nonconst();
  atm.initialize();
}

/**
 * @brief Execute one time step.
 * @param dt Time step size in seconds
 */
void emulator_atm_run(const int dt) {
  auto &atm = get_atm_emulator_nonconst();
  atm.run(dt);
}

/**
 * @brief Finalize and clean up the atmosphere emulator.
 */
void emulator_atm_finalize() {
  auto &atm = get_atm_emulator_nonconst();
  atm.finalize();

  emulator::cleanup_emulator_context();
}

/**
 * @brief Get number of local columns on this rank.
 */
int emulator_atm_get_num_local_cols() {
  const auto &atm = get_atm_emulator();
  return atm.get_num_local_cols();
}

/**
 * @brief Get total number of global columns.
 */
int emulator_atm_get_num_global_cols() {
  const auto &atm = get_atm_emulator();
  return atm.get_num_global_cols();
}

/**
 * @brief Get grid size in x-direction.
 */
int emulator_atm_get_nx() {
  const auto &atm = get_atm_emulator();
  return atm.get_nx();
}

/**
 * @brief Get grid size in y-direction.
 */
int emulator_atm_get_ny() {
  const auto &atm = get_atm_emulator();
  return atm.get_ny();
}

/**
 * @brief Get global IDs for local columns.
 * @param gids Output array (must be pre-allocated)
 */
void emulator_atm_get_local_cols_gids(int *gids) {
  const auto &atm = get_atm_emulator();
  atm.get_local_col_gids(gids);
}

/**
 * @brief Get latitude/longitude for local columns.
 * @param lat Output latitude array [radians]
 * @param lon Output longitude array [radians]
 */
void emulator_atm_get_cols_latlon(double *lat, double *lon) {
  const auto &atm = get_atm_emulator();
  atm.get_cols_latlon(lat, lon);
}

/**
 * @brief Get cell areas for local columns.
 * @param area Output area array
 */
void emulator_atm_get_cols_area(double *area) {
  const auto &atm = get_atm_emulator();
  atm.get_cols_area(area);
}

} // extern "C"
