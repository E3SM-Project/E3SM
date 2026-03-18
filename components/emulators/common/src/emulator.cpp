/**
 * @file emulator.cpp
 * @brief Implementation of the Emulator base class.
 */

#include "emulator.hpp"
#include <stdexcept>
#include <iostream>


namespace emulator {

Emulator::Emulator(EmulatorType type, int id, const std::string &name)
    : m_type(type), m_id(id), m_name(name), m_initialized(false),
      m_step_count(0) {}

void Emulator::initialize() {
  if (m_initialized) {
    throw std::runtime_error("Emulator already initialized");
  }
  init_impl();
  m_initialized = true;
}

void Emulator::run(int dt) {
  if (!m_initialized) {
    throw std::runtime_error("Emulator::run() called before initialize()");
  }
  run_impl(dt);
  m_step_count++;
}

void Emulator::finalize() {
  if (!m_initialized) {
    return; // Already finalized or never initialized
  }
  final_impl();
  m_initialized = false;
}


static std::string to_string(EmulatorType t) {
  switch (t) {
  case EmulatorType::ATM_COMP: return "ATM_COMP";
  case EmulatorType::OCN_COMP: return "OCN_COMP";
  case EmulatorType::ICE_COMP: return "ICE_COMP";
  case EmulatorType::LND_COMP: return "LND_COMP";
  default:                     return "UNKNOWN";
  }
}

void Emulator::print_info(std::ostream& os) const {
  os << "Emulator '" << m_name << "'\n";
  os << "  type          : " << to_string(m_type) << "\n";
  os << "  id            : " << m_id << "\n";
  os << "  initialized   : " << std::boolalpha << m_initialized << "\n";
  os << "  step_count    : " << m_step_count << "\n";

  // Grid summary via virtual getters
  int nx    = get_nx();
  int ny    = get_ny();
  int nloc  = get_num_local_cols();
  int nglob = get_num_global_cols();

  os << "  grid          : nx=" << nx
     << " ny=" << ny
     << " num_local_cols=" << nloc
     << " num_global_cols=" << nglob << "\n";

  // Hook for derived classes to print config / component-specific info
  print_extra_info(os);
}

} // namespace emulator
