/**
 * @file atm_f2c.cpp
 * @brief ATM-specific emulator_create implementation.
 *
 * Provides emulator_create for the ATM component so that the Fortran
 * atm_comp_mct module does not need to link against emulator_driver.
 * All other C API calls (emulator_init, emulator_run, etc.) are satisfied
 * by emulator_common/src/emulator_c_api.cpp.
 */

#include "atm.hpp"
#include "emulator_c_api.hpp"

#include <mpi.h>
#include <string>

extern "C" {

void *emulator_create(const char *kind, const EmulatorCreateConfig *cfg) {
  // This file only handles ATM; kind must be "atm".
  auto *atm = new emulator::EmulatorAtm();

  atm->create_instance(
      cfg->f_comm,
      cfg->comp_id,
      cfg->input_file ? cfg->input_file : "",
      cfg->log_file   ? cfg->log_file   : "",
      cfg->run_type,
      cfg->start_ymd,
      cfg->start_tod);

  return static_cast<void *>(atm);
}

} // extern "C"
