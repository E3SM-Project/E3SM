/**
 * @file atm_f2c.cpp
 * @brief C interface for the atmosphere emulator (Fortran-callable).
 *
 * Provides extern "C" functions callable from Fortran via ISO_C_BINDING.
 * These wrap the EmulatorAtm C++ class methods.
 *
 * Lifecycle:
 * 1. emulator_atm_create_instance
 * 2. emulator_atm_set_grid_data
 * 3. emulator_atm_init_coupling_indices
 * 4. emulator_atm_setup_coupling
 * 5. emulator_atm_init
 * 6. emulator_atm_run
 * 7. emulator_atm_finalize
 */

#include "atm.hpp"
#include "emulator_registry.hpp"
#include <fstream>
#include <string>
#include <mpi.h>

namespace {

static const std::string ATM_INSTANCE_NAME = "eatm_default";

// Minimal log stream – opened once on master rank, points to atm.log
static std::ofstream s_log_stream;
static bool s_log_open = false;

emulator::EmulatorAtm &get_atm() {
  return emulator::EmulatorRegistry::instance()
      .get_mut<emulator::EmulatorAtm>(ATM_INSTANCE_NAME);
}

const emulator::EmulatorAtm &get_atm_const() {
  return emulator::EmulatorRegistry::instance()
      .get<emulator::EmulatorAtm>(ATM_INSTANCE_NAME);
}

} // anonymous namespace

extern "C" {

void emulator_atm_create_instance(const int f_comm, const int comp_id,
                                  const char *input_file,
                                  const char *log_file,
                                  const int run_type,
                                  const int start_ymd,
                                  const int start_tod) {
  // Open log file on master rank (Fortran already opened it,
  // so we append).  Non-master ranks skip file logging.
  int rank = 0;
  MPI_Comm_rank(MPI_Comm_f2c(f_comm), &rank);
  if (rank == 0 && log_file && log_file[0] != '\0') {
    s_log_stream.open(log_file, std::ios::app);
    s_log_open = s_log_stream.is_open();
    if (s_log_open) {
      s_log_stream << "(eatm_f2c) create_instance: "
                   << "comp_id=" << comp_id
                   << " run_type=" << run_type
                   << " start_ymd=" << start_ymd
                   << " start_tod=" << start_tod
                   << std::endl;
    }
  }

  auto &reg = emulator::EmulatorRegistry::instance();
  auto &atm = reg.create<emulator::EmulatorAtm>(ATM_INSTANCE_NAME);

  std::string inf = input_file ? input_file : "";
  atm.create_instance(f_comm, comp_id, inf, run_type,
                      start_ymd, start_tod);

  if (s_log_open) {
    s_log_stream << "(eatm_f2c) create_instance done: "
                 << "nx=" << atm.get_nx()
                 << " ny=" << atm.get_ny()
                 << " ncols=" << atm.get_num_global_cols()
                 << std::endl;
  }
}

void emulator_atm_set_grid_data(const int nx, const int ny,
                                const int num_local_cols,
                                const int num_global_cols,
                                const int *col_gids,
                                const double *lat,
                                const double *lon,
                                const double *area) {
  get_atm().set_grid_data(nx, ny, num_local_cols, num_global_cols,
                          col_gids, lat, lon, area);
}

void emulator_atm_init_coupling_indices(const char *export_fields,
                                        const char *import_fields) {
  std::string exp_str = export_fields ? export_fields : "";
  std::string imp_str = import_fields ? import_fields : "";
  get_atm().init_coupling_indices(exp_str, imp_str);
}

void emulator_atm_setup_coupling(double *import_data,
                                 double *export_data,
                                 const int num_imports,
                                 const int num_exports,
                                 const int field_size) {
  get_atm().setup_coupling(import_data, export_data,
                           num_imports, num_exports, field_size);
}

void emulator_atm_init() {
  if (s_log_open) {
    s_log_stream << "(eatm_f2c) initialize starting"
                 << std::endl;
  }
  get_atm().initialize();
  if (s_log_open) {
    s_log_stream << "(eatm_f2c) initialize complete"
                 << std::endl;
  }
}

void emulator_atm_run(const int dt) {
  auto &atm = get_atm();
  // Log first few steps to confirm run loop is working
  if (s_log_open && atm.step_count() < 3) {
    s_log_stream << "(eatm_f2c) run step="
                 << atm.step_count()
                 << " dt=" << dt << std::endl;
  }
  atm.run(dt);
}

void emulator_atm_finalize() {
  if (s_log_open) {
    s_log_stream << "(eatm_f2c) finalize starting"
                 << std::endl;
  }
  get_atm().finalize();
  emulator::cleanup_emulator_registry();
  if (s_log_open) {
    s_log_stream << "(eatm_f2c) finalize complete"
                 << std::endl;
    s_log_stream.close();
    s_log_open = false;
  }
}

int emulator_atm_get_num_local_cols() {
  return get_atm_const().get_num_local_cols();
}

int emulator_atm_get_num_global_cols() {
  return get_atm_const().get_num_global_cols();
}

int emulator_atm_get_nx() {
  return get_atm_const().get_nx();
}

int emulator_atm_get_ny() {
  return get_atm_const().get_ny();
}

void emulator_atm_get_local_cols_gids(int *gids) {
  get_atm_const().get_local_col_gids(gids);
}

void emulator_atm_get_cols_latlon(double *lat, double *lon) {
  get_atm_const().get_cols_latlon(lat, lon);
}

void emulator_atm_get_cols_area(double *area) {
  get_atm_const().get_cols_area(area);
}

} // extern "C"
