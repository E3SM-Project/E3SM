/**
 * @file emulator.hpp
 * @brief Abstract base class for all E3SM emulators.
 */

#ifndef E3SM_EMULATOR_HPP
#define E3SM_EMULATOR_HPP

#include <string>

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
