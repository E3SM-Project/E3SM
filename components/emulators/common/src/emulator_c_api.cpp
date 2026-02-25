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

} // extern "C"
