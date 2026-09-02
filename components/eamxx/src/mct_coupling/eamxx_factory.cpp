/**
 * @file atm_factory.cpp
 * @brief Component-local emulator_create for the full E3SM build.
 *
 * In the integrated E3SM build, atm_comp_mct.F90 calls emulator_create
 * via the C interop interface declared in emulator_f2c_api.F90.
 * This file provides that symbol directly inside libemulatoratm so that
 * there is no circular static-library dependency on emulator_driver.
 *
 * In standalone/test builds, emulator_create is provided by
 * emulator_driver/src/emulator_factory.cpp which additionally supports
 * the full EmulatorRegistry and multi-component dispatch.
 *
 * This file is compiled only when NOT in STANDALONE_MODE (see
 * emulatoratm/CMakeLists.txt).
 */

#include "eamxx_atm.hpp"
#include <emulator_c_api.hpp>

#include <cstring>

extern "C" {

/**
 * @brief Create an EAMxx ATM instance.
 *
 * @param kind  Must be "atm".  Returns nullptr for any other kind.
 * @param cfg   Creation configuration.
 * @return Opaque pointer to the new EmulatorAtm, or nullptr on error.
 */
void *emulator_create(const char *kind, const EmulatorCreateConfig *cfg) {
  if (!kind || std::strcmp(kind, "atm") != 0) {
    return nullptr;
  }

  scream::EAMxxAtm::create_instance(cfg->comm, cfg->comp_id, cfg->input_file, cfg->log_file,
                                    cfg->run_type, cfg->start_ymd, cfg->start_tod,
                                    cfg->case_start_ymd, cfg->case_start_tod, cfg->calendar);

  return static_cast<void *>(atm);
}

} // extern "C"
