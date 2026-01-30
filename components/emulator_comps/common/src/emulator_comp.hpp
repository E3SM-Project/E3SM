/**
 * @file emulator_comp.hpp
 * @brief Abstract base class for all emulated E3SM components.
 */

#ifndef EMULATOR_COMP_HPP
#define EMULATOR_COMP_HPP

#include <string>

namespace emulator {

/**
 * @brief Enumeration of component types in E3SM.
 */
enum class CompType {
  ATM = 0, ///< Atmosphere component
  OCN = 1, ///< Ocean component
  ICE = 2, ///< Sea ice component
  LND = 3  ///< Land component
};

/**
 * @brief Abstract base class for all emulated E3SM components.
 *
 * Provides the common infrastructure for emulator components.
 * Derived classes implement the pure virtual methods for
 * component-specific behavior.
 */
class EmulatorComp {
public:
  explicit EmulatorComp(CompType type);
  virtual ~EmulatorComp() = default;

  // Lifecycle methods
  void create_instance(int comp_id, const std::string &name);
  void initialize();
  void run(int dt);
  void finalize();

  // Accessors
  CompType type() const { return m_type; }
  int comp_id() const { return m_comp_id; }
  const std::string &name() const { return m_name; }
  bool is_initialized() const { return m_initialized; }
  int step_count() const { return m_step_count; }

protected:
  // Virtual methods for derived classes
  virtual void init_impl() = 0;
  virtual void run_impl(int dt) = 0;
  virtual void final_impl() = 0;

  CompType m_type;
  int m_comp_id = -1;
  std::string m_name;
  bool m_initialized = false;
  int m_step_count = 0;
};

} // namespace emulator

#endif // EMULATOR_COMP_HPP
