/**
 * @file atm.hpp
 * @brief Atmosphere emulator component declaration.
 *
 * Defines the EmulatorAtm class which implements an AI-based atmosphere
 * component for E3SM. Inherits from the Emulator base class and adds
 * atmosphere-specific coupling, field management, and inference.
 */

#ifndef EMULATORATM_HPP
#define EMULATORATM_HPP

#include "emulator_c_api.hpp"
#include <emulator.hpp>
#include <memory>
#include <string>
#include <vector>

namespace emulator {

/**
 * @brief Atmosphere emulator component.
 *
 * Derived from Emulator, provides atmosphere-specific functionality:
 * - Coupling field mappings (x2a inputs, a2x outputs)
 * - AI model integration via configurable inference backends
 * - MCT interface for CIME integration
 *
 * Currently assumes a structured lat-lon grid. Grid dimensions (nx, ny)
 * are read from atm_in and used to compute the total global column count
 * as nx * ny. Lat/lon coordinates are stored and passed to MCT in degrees.
 */
class EmulatorAtm : public Emulator {
public:
  EmulatorAtm();
  ~EmulatorAtm() = default;

  void init_data() override;
protected:
  // Virtual methods from Emulator base
  void run_impl(int dt) override;
  void final_impl() override;
  void print_extra_info(std::ostream& os) const override {};

};

} // namespace emulator

#endif // EMULATORATM_HPP
