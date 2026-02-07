/**
 * @file atm.hpp
 * @brief Concrete atmosphere emulator for MCT coupling.
 *
 * Inherits from emulator::Emulator. Grid-agnostic: decomposition
 * is set externally by the F90 MCT layer, not computed internally.
 */

#ifndef EMULATOR_ATM_HPP
#define EMULATOR_ATM_HPP

#include "emulator.hpp"

namespace emulator {
namespace atm {

/**
 * @brief Atmosphere emulator component.
 *
 * Grid-agnostic design: local decomposition (lsize, gindex) is
 * set via set_decomposition() called from F90 MCT layer.
 * For standalone testing, use set_dummy_decomposition().
 */
class Atm : public Emulator {
public:
  /**
   * @brief Construct a new Atm emulator.
   *
   * @param fcomm Fortran MPI communicator handle (integer)
   * @param compid Component id assigned by the driver
   */
  Atm(int fcomm, int compid);

  /**
   * @brief Set decomposition from MCT layer.
   *
   * Called by F90 after obtaining decomposition from E3SM.
   *
   * @param lsize Local grid size on this PE
   * @param gindex Array of 1-based global indices (length lsize)
   */
  void set_decomposition(int lsize, const int *gindex);

  /**
   * @brief Set a simple dummy decomposition for standalone testing.
   *
   * Creates a contiguous 1-D decomposition across MPI ranks.
   *
   * @param global_size Total number of grid points globally
   */
  void set_dummy_decomposition(int global_size);

protected:
  // Emulator virtual overrides
  void init_impl() override;
  void run_impl(int dt) override;
  void final_impl() override;
};

} // namespace atm
} // namespace emulator

#endif // EMULATOR_ATM_HPP
