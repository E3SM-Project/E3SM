#include "emulator.hpp"
#include "emulator_c_api.hpp"
#include "emulator_registry.hpp"

#include "atm.hpp"
// #include "ocn.hpp"
// #include "lnd.hpp"

#include <cstring>
#include <fstream>
#include <mpi.h>
#include <string>

namespace {

// Minimal log stream – opened once on master rank, points to atm.log
static std::ofstream s_log_stream;
static bool s_log_open = false;


} // anonymous namespace
extern "C" {
void *emulator_create(const char *kind, const EmulatorCreateConfig *cfg) {
  auto &reg = emulator::EmulatorRegistry::instance();

  int rank = 0;
  MPI_Comm_rank(MPI_Comm_f2c(cfg->f_comm), &rank);
  if (rank == 0 && cfg->log_file && cfg->log_file[0] != '\0') {
    s_log_stream.open(cfg->log_file, std::ios::app);
    s_log_open = s_log_stream.is_open();
    if (s_log_open) {
      s_log_stream << "(emulator_create) create_instance: "
                   << "comp_id=" << cfg->comp_id
                   << " run_type=" << cfg->run_type
                   << " start_ymd=" << cfg->start_ymd
                   << " start_tod=" << cfg->start_tod << std::endl;
    }
  }

  if (std::strcmp(kind, "atm") == 0) {
    auto &atm = reg.create<emulator::EmulatorAtm>("eatm_" +
                                                  std::to_string(cfg->comp_id));
    // Configure it
    atm.create_instance(cfg->f_comm, cfg->comp_id,
                        cfg->input_file ? cfg->input_file : "", 
                        cfg->log_file ? cfg->log_file : "",
                        cfg->run_type,
                        cfg->start_ymd, cfg->start_tod);

    emulator::Emulator *base = &atm;
    if (s_log_open) {
      s_log_stream << "(emulator_factory) create_instance done!"<< std::endl;
    }
    // return opaque pointer to Fortran
    return static_cast<void *>(base);
  }

  // TODO: ocean, land…

  return nullptr;
}
} // extern "C"
