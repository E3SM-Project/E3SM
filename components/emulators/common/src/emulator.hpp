/**
 * @file emulator.hpp
 * @brief Abstract base class for all E3SM emulators.
 */

#ifndef E3SM_EMULATOR_HPP
#define E3SM_EMULATOR_HPP

#include <string>
#include "emulator_c_api.hpp"

namespace emulator {

/**
 * @brief Enumeration of emulator types in E3SM.
 */
enum class EmulatorType {
  ATM_COMP = 0, ///< Atmosphere component emulator
  OCN_COMP = 1, ///< Ocean component emulator
  ICE_COMP = 2, ///< Sea ice component emulator
  LND_COMP = 3  ///< Land component emulator
};

/**
 * @brief Abstract base class for all E3SM emulators.
 *
 * Provides the common infrastructure for emulators.
 * Derived classes implement the pure virtual methods for
 * emulator-specific behavior.
 */
class Emulator {
public:
  /**
   * @brief Construct a new Emulator.
   *
   * @param type Emulator type
   * @param id Emulator ID (-1 if unassigned)
   * @param name Emulator name (empty if unassigned)
   */
  explicit Emulator(EmulatorType type, int id = -1,
                    const std::string &name = "");
  virtual ~Emulator() = default;

  // Lifecycle methods
  void initialize();
  void run(int dt);
  void finalize();

  // Accessors
  EmulatorType type() const { return m_type; }
  int id() const { return m_id; }
  const std::string &name() const { return m_name; }
  bool is_initialized() const { return m_initialized; }
  int step_count() const { return m_step_count; }

  // New virtuals for grid / coupling
  virtual void set_grid_data(const EmulatorGridDesc& grid) = 0;
  virtual void setup_coupling(const EmulatorCouplingDesc& cpl) = 0;
  virtual void init_coupling_indices(const std::string &export_fields,
                             const std::string &import_fields) = 0;

  // Optionally: virtual accessors if Fortran needs them
  virtual int get_num_local_cols() const = 0;
  virtual int get_num_global_cols() const = 0;
  virtual int get_nx() const = 0;
  virtual int get_ny() const = 0;
  virtual void get_local_col_gids(int* gids) const = 0;
  virtual void get_cols_latlon(double* lat, double* lon) const = 0;
  virtual void get_cols_area(double* area) const = 0;

protected:
  // Virtual methods for derived classes
  virtual void init_impl() = 0;
  virtual void run_impl(int dt) = 0;
  virtual void final_impl() = 0;

  EmulatorType m_type;
  int m_id;
  std::string m_name;
  bool m_initialized = false;
  int m_step_count = 0;
};

} // namespace emulator

#endif // E3SM_EMULATOR_HPP
