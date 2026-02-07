/**
 * @file eocn_binding.cpp
 * @brief extern "C" functions bridging Fortran MCT calls to
 *        the C++ Ocn emulator via the EmulatorRegistry.
 *
 * Grid-agnostic: decomposition is passed in from F90, not computed here.
 */

#include "emulator_registry.hpp"
#include "ocn.hpp"

#include <cstring>

// Registry key for the singleton OCN emulator
static const char *EOCN_KEY = "eocn";

extern "C" {

/**
 * Create the Ocn emulator in the global registry.
 */
void eocn_create_c(int fcomm, int compid) {
  auto &reg = emulator::EmulatorRegistry::instance();
  reg.create<emulator::ocn::Ocn>(EOCN_KEY, fcomm, compid);
}

/**
 * Set the decomposition from F90 MCT layer.
 * gindex is an array of 1-based global indices.
 */
void eocn_set_decomposition_c(int lsize, const int *gindex) {
  auto &reg = emulator::EmulatorRegistry::instance();
  auto &ocn = reg.get_mut<emulator::ocn::Ocn>(EOCN_KEY);
  ocn.set_decomposition(lsize, gindex);
}

/**
 * Return the local grid size on this PE.
 */
int eocn_get_lsize_c() {
  auto &reg = emulator::EmulatorRegistry::instance();
  const auto &ocn = reg.get<emulator::ocn::Ocn>(EOCN_KEY);
  return ocn.lsize();
}

/**
 * Copy global indices into caller-supplied Fortran array.
 * Array must be pre-allocated to at least lsize elements.
 */
void eocn_get_gindex_c(int *gindex) {
  auto &reg = emulator::EmulatorRegistry::instance();
  const auto &ocn = reg.get<emulator::ocn::Ocn>(EOCN_KEY);
  std::memcpy(gindex, ocn.gindex().data(), ocn.lsize() * sizeof(int));
}

/**
 * Delegate to Emulator::initialize().
 */
void eocn_init_c() {
  auto &reg = emulator::EmulatorRegistry::instance();
  auto &ocn = reg.get_mut<emulator::ocn::Ocn>(EOCN_KEY);
  ocn.initialize();
}

/**
 * Delegate to Emulator::run(dt).
 */
void eocn_run_c(int dt) {
  auto &reg = emulator::EmulatorRegistry::instance();
  auto &ocn = reg.get_mut<emulator::ocn::Ocn>(EOCN_KEY);
  ocn.run(dt);
}

/**
 * Delegate to Emulator::finalize() and clean up registry.
 */
void eocn_final_c() {
  auto &reg = emulator::EmulatorRegistry::instance();
  if (reg.has(EOCN_KEY)) {
    auto &ocn = reg.get_mut<emulator::ocn::Ocn>(EOCN_KEY);
    ocn.finalize();
  }
  emulator::cleanup_emulator_registry();
}

} // extern "C"
