/**
 * @file eatm_binding.cpp
 * @brief extern "C" functions bridging Fortran MCT calls to
 *        the C++ AtmEmulator via the EmulatorRegistry.
 */

#include "eatm_atm_emulator.hpp"
#include "emulator_registry.hpp"

#include <cstring>

// Registry key for the singleton ATM emulator
static const char* EATM_KEY = "eatm";

extern "C" {

/**
 * Create the AtmEmulator in the global registry and compute
 * the grid decomposition.
 */
void eatm_create_c(int fcomm, int compid, int nxg, int nyg) {
  auto& reg = emulator::EmulatorRegistry::instance();
  auto& atm = reg.create<eatm::AtmEmulator>(
      EATM_KEY, fcomm, compid, nxg, nyg);
  atm.compute_decomposition();
}

/**
 * Return the local grid size on this PE.
 */
int eatm_get_lsize_c() {
  auto& reg = emulator::EmulatorRegistry::instance();
  const auto& atm = reg.get<eatm::AtmEmulator>(EATM_KEY);
  return atm.lsize();
}

/**
 * Copy grid data into caller-supplied Fortran arrays.
 * Arrays must be pre-allocated to at least lsize elements.
 */
void eatm_get_grid_data_c(int* gindex, double* lons, double* lats,
                           double* areas, double* masks, double* fracs) {
  auto& reg = emulator::EmulatorRegistry::instance();
  const auto& atm = reg.get<eatm::AtmEmulator>(EATM_KEY);
  int ls = atm.lsize();

  std::memcpy(gindex, atm.gindex().data(), ls * sizeof(int));
  std::memcpy(lons,   atm.lons().data(),   ls * sizeof(double));
  std::memcpy(lats,   atm.lats().data(),   ls * sizeof(double));
  std::memcpy(areas,  atm.areas().data(),  ls * sizeof(double));
  std::memcpy(masks,  atm.masks().data(),  ls * sizeof(double));
  std::memcpy(fracs,  atm.fracs().data(),  ls * sizeof(double));
}

/**
 * Delegate to Emulator::initialize().
 */
void eatm_init_c() {
  auto& reg = emulator::EmulatorRegistry::instance();
  auto& atm = reg.get_mut<eatm::AtmEmulator>(EATM_KEY);
  atm.initialize();
}

/**
 * Delegate to Emulator::run(dt).
 */
void eatm_run_c(int dt) {
  auto& reg = emulator::EmulatorRegistry::instance();
  auto& atm = reg.get_mut<eatm::AtmEmulator>(EATM_KEY);
  atm.run(dt);
}

/**
 * Delegate to Emulator::finalize() and clean up registry.
 */
void eatm_final_c() {
  auto& reg = emulator::EmulatorRegistry::instance();
  if (reg.has(EATM_KEY)) {
    auto& atm = reg.get_mut<eatm::AtmEmulator>(EATM_KEY);
    atm.finalize();
  }
  emulator::cleanup_emulator_registry();
}

} // extern "C"
