

# File emulator\_context.hpp

[**File List**](files.md) **>** [**common**](dir_03a82009d675450cd7b26f7887b718e6.md) **>** [**src**](dir_cf65ed25a3cb4ff7a3950496111c5548.md) **>** [**emulator\_context.hpp**](emulator__context_8hpp.md)

[Go to the documentation of this file](emulator__context_8hpp.md)


```C++


#ifndef EMULATOR_CONTEXT_HPP
#define EMULATOR_CONTEXT_HPP

#include <any>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeindex>

namespace emulator {

class EmulatorContext {
public:
  static EmulatorContext &singleton() {
    static EmulatorContext c;
    return c;
  }

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

  template <typename T> bool has() const {
    auto key = getKey<T>();
    return m_objects.find(key) != m_objects.end();
  }

  void clean_up() { m_objects.clear(); }

private:
  using key_type = std::type_index;

  EmulatorContext() = default;

  template <typename T> static key_type getKey() {
    return std::type_index(typeid(T));
  }

  std::map<key_type, std::any> m_objects; 
};

inline void cleanup_emulator_context() {
  EmulatorContext::singleton().clean_up();
}

} // namespace emulator

#endif // EMULATOR_CONTEXT_HPP
```


