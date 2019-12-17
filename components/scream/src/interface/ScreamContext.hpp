#ifndef SCREAM_CONTEXT_HPP
#define SCREAM_CONTEXT_HPP

#include "share/util/scream_std_any.hpp"
#include <map>
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
    scream_require_msg(m_objects.find(key)==m_objects.end(),
                       "Error! Object with key '" + key.name() + "' was already created in the scream context.\n");

    auto it = m_objects.emplace(key,args...);

    scream_require_msg(it.second, "Error! Something went wrong while creating object with key '" + key.name() + "'.\n");
    return *it.first;

  }

  template<typename T>
  const T& get () const {
    auto key = getKey<T>();
    scream_require_msg(m_objects.find(key)!=m_objects.end(),
                       "Error! Object with key '" + key.name() + "' not found in the scream context.\n");
    return m_objects.at(key);
  }

  template<typename T>
  T& getNonConst () const {
    auto key = getKey<T>();
    scream_require_msg(m_objects.find(key)!=m_objects.end(),
                       "Error! Object with key '" + key.name() + "' not found in the scream context.\n");
    return m_objects.at(key);
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

  std::map<key_type,util::any>  m_objects;
};

} // namespace scream

#endif // SCREAM_CONTEXT_HPP
