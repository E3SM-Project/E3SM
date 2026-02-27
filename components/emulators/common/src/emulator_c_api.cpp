#include "emulator_c_api.hpp"
#include "emulator.hpp"

#include <cstring>
#include <iostream>
extern "C" {

void emulator_set_grid_data(void* handle,
                            const EmulatorGridDesc* grid) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->set_grid_data(*grid);
}

void emulator_setup_coupling(void* handle,
                             EmulatorCouplingDesc* cpl) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->setup_coupling(*cpl);
}

void emulator_init(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->initialize();
}

void emulator_run(void* handle, int dt) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->run(dt);
}

void emulator_finalize(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->finalize();
}

void emulator_print_info(void *handle){
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->print_info(std::cout);
}
void emulator_init_coupling_indices(void* handle, const char* import_fields, const char* export_fields){
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->init_coupling_indices(std::string(export_fields), std::string(import_fields));
}

int emulator_get_num_local_cols(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  return emu->get_num_local_cols();
}

int emulator_get_num_global_cols(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  return emu->get_num_global_cols();
}

int emulator_get_nx(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  return emu->get_nx();
}

int emulator_get_ny(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  return emu->get_ny();
}

void emulator_get_local_col_gids(void* handle, int* gids) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->get_local_col_gids(gids);
}

void emulator_get_cols_latlon(void* handle, double* lat, double* lon) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->get_cols_latlon(lat, lon);
}

void emulator_get_cols_area(void* handle, double* area) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->get_cols_area(area);
}

} // extern "C"
