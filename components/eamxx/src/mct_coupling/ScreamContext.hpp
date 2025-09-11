#ifndef SCREAM_CONTEXT_HPP
#define SCREAM_CONTEXT_HPP

#include <map>
#include <any>
#include <typeindex>

namespace scream
{

class ScreamContext
{
public:

  static ScreamContext& singleton() {
    static ScreamContext c;
    return c;
  }

  template<typename T,typename... Args>
  T& create (Args... args) {
    auto key = getKey<T>();
    EKAT_REQUIRE_MSG(m_objects.find(key)==m_objects.end(),"Error! Object with key '" + (std::string)key.name() + "' was already created in the scream context.\n");

    auto& obj = m_objects[key];
    obj = std::make_any<T>(args...);

    return std::any_cast<T&>(obj);
  }

  template<typename T>
  const T& get () const {
    auto key = getKey<T>();
    EKAT_REQUIRE_MSG(m_objects.find(key)!=m_objects.end(),
                       "Error! Object with key '" + (std::string)key.name() + "' not found in the scream context.\n");
    const auto& obj = m_objects.at(key);

    return std::any_cast<const T&>(obj);
  }

  template<typename T>
  T& getNonConst () {
    auto key = getKey<T>();
    EKAT_REQUIRE_MSG(m_objects.find(key)!=m_objects.end(),
                       "Error! Object with key '" + (std::string)key.name() + "' not found in the scream context.\n");
    auto& obj = m_objects.at(key);

    return std::any_cast<T&>(obj);
  }

  void clean_up () {
    m_objects.clear();
  }

private:

  using key_type = std::type_index;

  ScreamContext () = default;

  // Rely on typeid to get a unique name for each type
  template<typename T>
  static key_type getKey() { return std::type_index(typeid(T)); }

  std::map<key_type,std::any>  m_objects;
};

inline void cleanup_singleton() {
  ScreamContext::singleton().clean_up();
}

} // namespace scream

#endif // SCREAM_CONTEXT_HPP
