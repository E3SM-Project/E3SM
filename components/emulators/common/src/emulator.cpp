/**
 * @file emulator.cpp
 * @brief Implementation of the Emulator base class.
 */

#include "emulator.hpp"
#include <stdexcept>

namespace emulator {

Emulator::Emulator(EmulatorType type, int fcomm, int id,
                   const std::string &name)
    : m_type(type), m_id(id), m_name(name), m_initialized(false),
      m_step_count(0), m_lsize(0) {
  // Convert Fortran communicator to C communicator
  m_comm = MPI_Comm_f2c(fcomm);
  if (m_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(m_comm, &m_rank);
    MPI_Comm_size(m_comm, &m_nprocs);
  } else {
    m_rank = 0;
    m_nprocs = 1;
  }
}

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

} // namespace emulator
