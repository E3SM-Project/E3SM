/**
 * @file emulator_context.hpp
 * @brief Singleton context for managing emulator component instances.
 *
 * Inspired by EAMxx's ScreamContext, this provides a type-safe registry
 * for creating and retrieving emulator component instances.
 */

#ifndef EMULATOR_CONTEXT_HPP
#define EMULATOR_CONTEXT_HPP

#include <any>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeindex>

namespace emulator {

/**
 * @brief Singleton context for managing emulator component instances.
 *
 * Provides type-safe storage and retrieval of component instances,
 * allowing the Fortran interface to create components and later
 * retrieve them for operations.
 *
 * ## Usage
 * ```cpp
 * // Create a component
 * auto& atm = EmulatorContext::singleton().create<EmulatorAtm>();
 *
 * // Later, retrieve it
 * auto& atm = EmulatorContext::singleton().getNonConst<EmulatorAtm>();
 *
 * // Check if it exists
 * if (EmulatorContext::singleton().has<EmulatorAtm>()) { ... }
 *
 * // Clean up
 * cleanup_emulator_context();
 * ```
 *
 * @note Only one instance of each type can exist in the context.
 */
class EmulatorContext {
public:
  /**
   * @brief Get the singleton instance.
   * @return Reference to the global EmulatorContext
   */
  static EmulatorContext &singleton() {
    static EmulatorContext c;
    return c;
  }

  /**
   * @brief Create and register a new component instance.
   *
   * Creates an instance of type T with the given constructor arguments
   * and stores it in the context.
   *
   * @tparam T Component type to create
   * @tparam Args Constructor argument types
   * @param args Arguments to pass to T's constructor
   * @return Reference to the newly created instance
   * @throws std::runtime_error if an instance of T already exists
   */
  template <typename T, typename... Args> T &create(Args &&...args) {
    auto key = getKey<T>();
    if (m_objects.find(key) != m_objects.end()) {
      throw std::runtime_error(
          "Error! Object with key '" + std::string(key.name()) +
          "' was already created in the emulator context.\n");
    }

    auto ptr = std::make_shared<T>(std::forward<Args>(args)...);
    m_objects[key] = std::any(ptr);

    return *ptr;
  }

  /**
   * @brief Get a const reference to an existing component.
   *
   * @tparam T Component type to retrieve
   * @return Const reference to the component
   * @throws std::runtime_error if no instance of T exists
   */
  template <typename T> const T &get() const {
    auto key = getKey<T>();
    auto it = m_objects.find(key);
    if (it == m_objects.end()) {
      throw std::runtime_error("Error! Object with key '" +
                               std::string(key.name()) +
                               "' not found in the emulator context.\n");
    }
    return *std::any_cast<std::shared_ptr<T>>(it->second);
  }

  /**
   * @brief Get a non-const reference to an existing component.
   *
   * @tparam T Component type to retrieve
   * @return Reference to the component
   * @throws std::runtime_error if no instance of T exists
   */
  template <typename T> T &getNonConst() {
    auto key = getKey<T>();
    auto it = m_objects.find(key);
    if (it == m_objects.end()) {
      throw std::runtime_error("Error! Object with key '" +
                               std::string(key.name()) +
                               "' not found in the emulator context.\n");
    }
    return *std::any_cast<std::shared_ptr<T>>(it->second);
  }

  /**
   * @brief Check if a component of type T exists.
   *
   * @tparam T Component type to check
   * @return true if an instance of T exists in the context
   */
  template <typename T> bool has() const {
    auto key = getKey<T>();
    return m_objects.find(key) != m_objects.end();
  }

  /**
   * @brief Remove all objects from the context.
   *
   * Should be called during shutdown to release resources.
   */
  void clean_up() { m_objects.clear(); }

private:
  using key_type = std::type_index;

  EmulatorContext() = default;

  template <typename T> static key_type getKey() {
    return std::type_index(typeid(T));
  }

  std::map<key_type, std::any> m_objects; ///< Type-indexed object storage
};

/**
 * @brief Convenience function to clean up the global emulator context.
 *
 * Equivalent to EmulatorContext::singleton().clean_up().
 */
inline void cleanup_emulator_context() {
  EmulatorContext::singleton().clean_up();
}

} // namespace emulator

#endif // EMULATOR_CONTEXT_HPP
