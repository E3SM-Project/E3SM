/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_CONTEXT_HPP
#define HOMMEXX_CONTEXT_HPP

#include <string>
#include <map>
#include <memory>

#include "ErrorDefs.hpp"
#include "utilities/StdMeta.hpp"

namespace Homme {

/* A Context manages resources that are morally singletons. Context is
 * meant to have two roles. First, a Context singleton is the only singleton in
 * the program. Second, a context need not be a singleton, and each Context
 * object can have different Elements, ReferenceElement, etc., objects. (That
 * probably isn't needed, but Context immediately supports it.)
 *
 * Finally, Context has two singleton functions: singleton(), which returns
 * Context&, and finalize_singleton(). The second is called in a unit test exe
 * main before Kokkos::finalize().
 */
class Context {
public:

  // Getters for a managed object.
  template<typename ConcreteType>
  bool has () const;

  // Getters for a managed object.
  template<typename ConcreteType, typename... Args>
  ConcreteType& get (Args&&... args);

  template<typename ConcreteType, typename... Args>
  std::shared_ptr<ConcreteType> get_ptr (Args&&... args);

  // Exactly one singleton.
  static Context& singleton();

  static void finalize_singleton();
private:

  std::map<std::string,Homme::any> m_members;

  // Clear the objects Context manages.
  void clear();
};

// ==================== IMPLEMENTATION =================== //

template<typename ConcreteType>
bool Context::has () const {
  const std::string& name = typeid(ConcreteType).name();
  auto it = m_members.find(name);
  return it!=m_members.end();
}

template<typename ConcreteType, typename... Args>
ConcreteType& Context::get (Args&&... args) {
  return *get_ptr<ConcreteType>(args...);
}

template<typename ConcreteType, typename... Args>
std::shared_ptr<ConcreteType> Context::get_ptr(Args&&... args) {
  if (!has<ConcreteType>()) {
    // This is needed for emplacing a type whose constructor takes no arguments.
    // We could do emplace(name,ConcreteType()), but then we would be assuming
    // that ConcreteType *has* a move constructor. This implementation here is
    // probably the most cumbersome, but also the safest.
    auto it_bool = m_members.emplace(std::piecewise_construct,
                                     std::make_tuple(typeid(ConcreteType).name()),
                                     std::make_tuple());
    Errors::runtime_check(it_bool.second, "Error! Something went wrong when inserting a new element in the context. "
                                          "This is an internal error. Please, contact developers.\n", -1);
    it_bool.first->second.reset<ConcreteType>(args...);
  }

  const std::string& name = typeid(ConcreteType).name();
  auto it = m_members.find(name);
  Errors::runtime_check(it!=m_members.end(), "Error! Context member '" + name + "' not found. This is an internal bug. Please, contact developers.\n", -1);

  return any_ptr_cast<ConcreteType>(it->second);
}

} // namespace Homme

#endif // HOMMEXX_CONTEXT_HPP
