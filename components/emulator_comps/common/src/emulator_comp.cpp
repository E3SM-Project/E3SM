/**
 * @file emulator_comp.cpp
 * @brief Implementation of the EmulatorComp base class.
 */

#include "emulator_comp.hpp"
#include <stdexcept>

namespace emulator {

EmulatorComp::EmulatorComp(CompType type)
    : m_type(type), m_comp_id(-1), m_initialized(false), m_step_count(0) {}

void EmulatorComp::create_instance(int comp_id, const std::string &name) {
  m_comp_id = comp_id;
  m_name = name;
}

void EmulatorComp::initialize() {
  if (m_initialized) {
    throw std::runtime_error("EmulatorComp already initialized");
  }
  init_impl();
  m_initialized = true;
}

void EmulatorComp::run(int dt) {
  if (!m_initialized) {
    throw std::runtime_error("EmulatorComp::run() called before initialize()");
  }
  run_impl(dt);
  m_step_count++;
}

void EmulatorComp::finalize() {
  if (!m_initialized) {
    return; // Already finalized or never initialized
  }
  final_impl();
  m_initialized = false;
}

} // namespace emulator
