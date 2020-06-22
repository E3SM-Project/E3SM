#ifndef SCREAM_PARAMETER_LIST_HPP
#define SCREAM_PARAMETER_LIST_HPP

#include "util/scream_std_meta.hpp"
#include "scream_assert.hpp"

#include <map>

namespace scream {

class ParameterList {
public:

  ParameterList () = default;
  explicit ParameterList (const std::string& name) : m_name(name) {}

  ParameterList& operator= (const ParameterList&) = default;

  const std::string& name () const { return m_name; }

  template<typename T>
  T& get (const std::string& name);

  template<typename T>
  const T& get (const std::string& name) const;

  template<typename T>
  T& get (const std::string& name, const T& def_value);

  template<typename T>
  void set (const std::string& name, const T& value);

  ParameterList& sublist (const std::string& name);

  const ParameterList& sublist (const std::string& name) const { return m_sublists.at(name); }

  bool isParameter (const std::string& name) const { return m_params.find(name)!=m_params.end(); }
  bool isSublist   (const std::string& name) const { return m_sublists.find(name)!=m_sublists.end(); }

  void print (std::ostream& out = std::cout, const int indent = 0) const;
private:

  std::string                           m_name;
  std::map<std::string,util::any>       m_params;
  std::map<std::string,ParameterList>   m_sublists;
};

// ====================== IMPLEMENTATION ===================== //

template<typename T>
inline T& ParameterList::get (const std::string& name) {
  // Check entry exists
  error::runtime_check ( isParameter(name),
                        "Error! Key '" + name + "' not found in parameter list '" + m_name + "'.\n");

  return util::any_cast<T>(m_params[name]);
}

template<typename T>
inline const T& ParameterList::get (const std::string& name) const {
  error::runtime_check ( isParameter(name),
                        "Error! Key '" + name + "' not found in parameter list '" + m_name + "'.\n");

  return util::any_cast<T>(m_params.at(name));
}

template<typename T>
inline T& ParameterList::get (const std::string& name, const T& def_value) {
  if ( !isParameter(name) ) {
    m_params[name].template reset<T>(def_value);
  }
  return get<T>(name);
}

template<typename T>
inline void ParameterList::set (const std::string& name, const T& value) {
  if ( !isParameter(name) ) {
    m_params[name].template reset<T>(value);
  } else {
    get<T>(name) = value;
  }
}

} // namespace scream

#endif // SCREAM_PARAMETER_LIST_HPP
