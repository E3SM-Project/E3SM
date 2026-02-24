#include "emulator_c_api.hpp"
#include "emulator.hpp"

#include <cstring>
extern "C" {

void emulator_set_grid_data(void* handle,
                            const EmulatorGridDesc* grid) {
  auto* emu = static_cast<emulator::Emulator*>(handle);  // <-- cast #2
  emu->set_grid_data(*grid);  // virtual → EmulatorAtm::set_grid_data
}

void emulator_setup_coupling(void* handle,
                             EmulatorCouplingDesc* cpl) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->setup_coupling(*cpl);
}

void emulator_init(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->initialize();          // calls init_impl()
}

void emulator_run(void* handle, int dt) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->run(dt);               // calls run_impl()
}

void emulator_finalize(void* handle) {
  auto* emu = static_cast<emulator::Emulator*>(handle);
  emu->finalize();            // calls final_impl()
}


} // extern "C"
