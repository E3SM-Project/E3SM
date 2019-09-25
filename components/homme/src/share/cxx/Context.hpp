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
#include <functional>

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

  // Setters for a managed object.
  template<typename ConcreteType, typename... Args>
  ConcreteType& create (Args&&... args);

  template<typename ConcreteType>
  void create_ref (ConcreteType& src);

  // More relaxed than create, it won't throw if the
  // object already exists, and simply return it.
  // NOTE: this is allows more flexibility in the cxx-f90
  //       interfaces, in that you have more freedom in
  //       the order of certain function calls.
  //       With 'create', you have to either check if an
  //       object was already created (and create it if it wasn't),
  //       or make sure that the function that calls 'create<T>'
  //       is *always* invoked *before* those that call 'get<T>'.
  template<typename ConcreteType, typename... Args>
  ConcreteType& create_if_not_there (Args&&... args);

  // Getters for a managed object.
  template<typename ConcreteType>
  ConcreteType& get () const;

  template<typename ConcreteType>
  std::shared_ptr<ConcreteType> get_ptr () const;

  // Exactly one singleton.
  static Context& singleton();

  static void finalize_singleton();
private:

  std::map<std::string,Homme::any> m_members;
  std::map<std::string, bool>      m_is_ref_wrapper;

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
ConcreteType& Context::create_if_not_there (Args&&... args) {
  // This is needed for emplacing a type whose constructor takes no arguments.
  // We could do emplace(name,ConcreteType()), but then we would be assuming
  // that ConcreteType *has* a move constructor. This implementation here is
  // probably the most cumbersome, but also the safest.
  auto it_bool = m_members.emplace(typeid(ConcreteType).name(),Homme::any());
  if (it_bool.second) {
    // We created it, so init it.
    it_bool.first->second.reset<ConcreteType>(args...);
    m_is_ref_wrapper[typeid(ConcreteType).name()] = false;
  }
  return *any_ptr_cast<ConcreteType>(it_bool.first->second);
}

template<typename ConcreteType, typename... Args>
ConcreteType& Context::create (Args&&... args) {
  Errors::runtime_check(!has<ConcreteType>(),
                        "Error! An object for the concrete type " + std::string(typeid(ConcreteType).name()) +
                        " is already stored. The 'Context' class does not allow overwriting or duplicates.\n");
  // This is needed for emplacing a type whose constructor takes no arguments.
  // We could do emplace(name,ConcreteType()), but then we would be assuming
  // that ConcreteType *has* a move constructor. This implementation here is
  // probably the most cumbersome, but also the safest.
  auto it_bool = m_members.emplace(typeid(ConcreteType).name(),Homme::any());
  Errors::runtime_check(it_bool.second, "Error! Something went wrong when inserting a new element in the context. "
                                        "This is an internal error. Please, contact developers.\n", -1);
  it_bool.first->second.reset<ConcreteType>(args...);
  m_is_ref_wrapper[typeid(ConcreteType).name()] = false;

  return *any_ptr_cast<ConcreteType>(it_bool.first->second);
}

template<typename ConcreteType>
void Context::create_ref (ConcreteType& src) {
  // This is needed for emplacing a type whose constructor takes no arguments.
  // We could do emplace(name,ConcreteType()), but then we would be assuming
  // that ConcreteType *has* a move constructor. This implementation here is
  // probably the most cumbersome, but also the safest.
  auto it_bool = m_members.emplace(typeid(ConcreteType).name(),Homme::any());
  Errors::runtime_check(it_bool.second, "Error! Something went wrong when inserting a new element in the context. "
                                        "This is an internal error. Please, contact developers.\n", -1);

  auto ref = std::ref(src);
  it_bool.first->second.reset<std::reference_wrapper<ConcreteType>>(ref);
  m_is_ref_wrapper[typeid(ConcreteType).name()] = true;
}

template<typename ConcreteType>
ConcreteType& Context::get () const {
  const std::string& name = typeid(ConcreteType).name();
  auto it = m_members.find(name);
  Errors::runtime_check(it!=m_members.end(), "Error! Context member '" + name + "' not found.\n", -1);
  if (m_is_ref_wrapper.at(name)) {
    auto ref = any_ptr_cast<std::reference_wrapper<ConcreteType>>(it->second);
    return ref->get();
  } else {
    return *get_ptr<ConcreteType>();
  }
}

template<typename ConcreteType>
std::shared_ptr<ConcreteType> Context::get_ptr() const {
  const std::string& name = typeid(ConcreteType).name();
  auto it = m_members.find(name);
  Errors::runtime_check(it!=m_members.end(), "Error! Context member '" + name + "' not found.\n", -1);
  Errors::runtime_check(!m_is_ref_wrapper.at(name),
                        "Error! Context member '" + name + "' is only available as a reference.\n", -1);

  return any_ptr_cast<ConcreteType>(it->second);
}

} // namespace Homme

#endif // HOMMEXX_CONTEXT_HPP
