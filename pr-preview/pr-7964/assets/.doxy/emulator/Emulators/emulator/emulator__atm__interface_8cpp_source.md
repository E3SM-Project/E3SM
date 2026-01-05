

# File emulator\_atm\_interface.cpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**emulator\_atm\_interface.cpp**](emulator__atm__interface_8cpp.md)

[Go to the documentation of this file](emulator__atm__interface_8cpp.md)


```C++


#include "../../common/src/emulator_comp.hpp"
#include "../../common/src/emulator_context.hpp"
#include "emulator_atm.hpp"

namespace {

emulator::EmulatorAtm &get_atm_emulator_nonconst() {
  auto &ctx = emulator::EmulatorContext::singleton();
  return ctx.getNonConst<emulator::EmulatorAtm>();
}

const emulator::EmulatorAtm &get_atm_emulator() {
  const auto &ctx = emulator::EmulatorContext::singleton();
  return ctx.get<emulator::EmulatorAtm>();
}

} // anonymous namespace

extern "C" {

// ===========================================================================
// Atmosphere Emulator C Interface (callable from Fortran)
// ===========================================================================

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

void emulator_atm_set_grid_data(const int nx, const int ny,
                                const int num_local_cols,
                                const int num_global_cols, const int *col_gids,
                                const double *lat, const double *lon,
                                const double *area) {
  auto &atm = get_atm_emulator_nonconst();
  atm.set_grid_data(nx, ny, num_local_cols, num_global_cols, col_gids, lat, lon,
                    area);
}

void emulator_atm_init_coupling_indices(const char *export_fields,
                                        const char *import_fields) {
  auto &atm = get_atm_emulator_nonconst();
  std::string exp_str = export_fields ? export_fields : "";
  std::string imp_str = import_fields ? import_fields : "";
  atm.init_coupling_indices(exp_str, imp_str);
}

void emulator_atm_setup_coupling(double *import_data, double *export_data,
                                 const int num_imports, const int num_exports,
                                 const int field_size) {
  auto &atm = get_atm_emulator_nonconst();
  atm.setup_coupling(import_data, export_data, num_imports, num_exports,
                     field_size);
}

void emulator_atm_init() {
  auto &atm = get_atm_emulator_nonconst();
  atm.initialize();
}

void emulator_atm_run(const int dt) {
  auto &atm = get_atm_emulator_nonconst();
  atm.run(dt);
}

void emulator_atm_finalize() {
  auto &atm = get_atm_emulator_nonconst();
  atm.finalize();

  emulator::cleanup_emulator_context();
}

int emulator_atm_get_num_local_cols() {
  const auto &atm = get_atm_emulator();
  return atm.get_num_local_cols();
}

int emulator_atm_get_num_global_cols() {
  const auto &atm = get_atm_emulator();
  return atm.get_num_global_cols();
}

int emulator_atm_get_nx() {
  const auto &atm = get_atm_emulator();
  return atm.get_nx();
}

int emulator_atm_get_ny() {
  const auto &atm = get_atm_emulator();
  return atm.get_ny();
}

void emulator_atm_get_local_cols_gids(int *gids) {
  const auto &atm = get_atm_emulator();
  atm.get_local_col_gids(gids);
}

void emulator_atm_get_cols_latlon(double *lat, double *lon) {
  const auto &atm = get_atm_emulator();
  atm.get_cols_latlon(lat, lon);
}

void emulator_atm_get_cols_area(double *area) {
  const auto &atm = get_atm_emulator();
  atm.get_cols_area(area);
}

} // extern "C"
```


